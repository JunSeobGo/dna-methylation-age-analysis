#!/usr/bin/env Rscript

# 450K 전체 beta matrix를 메모리에 올리지 않고 DNAm age 입력용 21k CpG만 추출한다.

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

beta_path <- file.path(project_dir, "data", "raw", "GSE40279", "GSE40279_average_beta.txt.gz")
manifest_path <- file.path(project_dir, "data", "processed", "GSE40279_sample_manifest.csv")
annotation_path <- file.path(project_dir, "data", "raw", "AdditionalFile22probeAnnotation21kdatMethUsed.csv")
output_path <- file.path(project_dir, "data", "processed", "GSE40279_beta_21k.tsv")
summary_path <- file.path(project_dir, "outputs", "gse40279_21k_extraction_summary.csv")

for (path in c(beta_path, manifest_path, annotation_path)) {
  if (!file.exists(path)) stop("필수 입력 파일을 찾을 수 없습니다: ", path)
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
annotation <- read.csv(annotation_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"Name" %in% names(annotation)) stop("21k CpG annotation에 Name 열이 없습니다.")
if (!"beta_column_name" %in% names(manifest)) stop("표본 매니페스트에 beta_column_name 열이 없습니다.")

required_probes <- unique(as.character(annotation$Name))
required_probes <- required_probes[!is.na(required_probes) & nzchar(required_probes)]
if (anyDuplicated(required_probes)) stop("21k CpG annotation에 중복 probe ID가 있습니다.")

input_connection <- gzfile(beta_path, open = "rt")
on.exit(close(input_connection), add = TRUE)
header <- strsplit(readLines(input_connection, n = 1L), "\t", fixed = TRUE)[[1]]
if (length(header) != nrow(manifest) + 1L || header[[1]] != "ID_REF") {
  stop("beta matrix의 열 구조가 표본 매니페스트와 일치하지 않습니다.")
}
if (!identical(header[-1], manifest$beta_column_name)) {
  stop("beta matrix의 표본 열 순서가 표본 매니페스트와 일치하지 않습니다.")
}

message("입력 검증 완료: ", length(required_probes), "개 CpG, ", nrow(manifest), "개 표본")
if (dry_run) quit(status = 0)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
temporary_output <- paste0(output_path, ".part")
if (file.exists(temporary_output)) unlink(temporary_output)
output_connection <- file(temporary_output, open = "wt", encoding = "UTF-8")
writeLines(paste(header, collapse = "\t"), output_connection)

found_probes <- setNames(rep(FALSE, length(required_probes)), required_probes)
input_row_count <- 0L
selected_row_count <- 0L
non_numeric_count <- 0L
out_of_range_count <- 0L
beta_min <- Inf
beta_max <- -Inf

repeat {
  lines <- readLines(input_connection, n = 1000L)
  if (!length(lines)) break
  input_row_count <- input_row_count + length(lines)
  probe_ids <- sub("\t.*$", "", lines)
  keep <- probe_ids %in% required_probes
  if (!any(keep)) next

  selected_lines <- lines[keep]
  selected_ids <- probe_ids[keep]
  for (index in seq_along(selected_lines)) {
    fields <- strsplit(selected_lines[[index]], "\t", fixed = TRUE)[[1]]
    if (length(fields) != length(header)) stop("CpG 행의 열 수가 header와 다릅니다: ", selected_ids[[index]])
    beta_text <- fields[-1]
    beta_values <- suppressWarnings(as.numeric(beta_text))
    non_numeric_count <- non_numeric_count + sum(is.na(beta_values) & !is.na(beta_text) & nzchar(beta_text))
    finite_values <- beta_values[is.finite(beta_values)]
    if (length(finite_values)) {
      beta_min <- min(beta_min, min(finite_values))
      beta_max <- max(beta_max, max(finite_values))
      out_of_range_count <- out_of_range_count + sum(finite_values < 0 | finite_values > 1)
    }
  }

  writeLines(selected_lines, output_connection)
  found_probes[selected_ids] <- TRUE
  selected_row_count <- selected_row_count + length(selected_lines)
}

close(output_connection)
missing_probes <- names(found_probes)[!found_probes]
if (length(missing_probes)) {
  unlink(temporary_output)
  stop("beta matrix에 없는 필수 CpG가 있습니다: ", length(missing_probes), "개")
}
if (selected_row_count != length(required_probes)) {
  unlink(temporary_output)
  stop("추출 CpG 행 수가 annotation과 일치하지 않습니다.")
}
if (non_numeric_count > 0L || out_of_range_count > 0L) {
  unlink(temporary_output)
  stop("beta 값 검증 실패: non_numeric=", non_numeric_count, ", out_of_range=", out_of_range_count)
}

if (file.exists(output_path)) unlink(output_path)
if (!file.rename(temporary_output, output_path)) stop("추출 파일 이름 변경에 실패했습니다.")

summary_table <- data.frame(
  dataset_id = "GSE40279",
  input_rows = input_row_count,
  selected_probes = selected_row_count,
  samples = nrow(manifest),
  beta_min = beta_min,
  beta_max = beta_max,
  non_numeric_count = non_numeric_count,
  out_of_range_count = out_of_range_count,
  stringsAsFactors = FALSE
)
write.csv(summary_table, summary_path, row.names = FALSE)
message("21k CpG 추출 완료: ", output_path)
print(summary_table)
