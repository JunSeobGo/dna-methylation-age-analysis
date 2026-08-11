#!/usr/bin/env Rscript

# SATSA 공개 beta 5개 파일에서 잠근 447명과 869개 모델 CpG만 행 단위로 추출한다.
# 성능을 계산하지 않으며, 사전 품질 게이트를 통과한 경우에만 외부 검증 bundle을 저장한다.

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
download <- "--download" %in% args
extract <- "--extract" %in% args
confirm_extraction <- "--confirm-extraction" %in% args
force <- "--force" %in% args
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

accession <- "S-BSST1206"
expected_source_rows <- 250816L
expected_source_samples <- 1469L
expected_selected_samples <- 447L
expected_model_features <- 869L
minimum_feature_coverage <- 0.95
maximum_sample_missing_rate <- 0.05
chunk_lines <- 1000L

download_manifest_path <- file.path(project_dir, "config", "dataset_download_manifest.csv")
sample_manifest_path <- file.path(project_dir, "data", "processed", "satsa_external_sample_manifest.csv")
candidate_model_path <- file.path(project_dir, "data", "processed", "locked_age_density_candidate.rds")
candidate_lock_path <- file.path(project_dir, "outputs", "age_density_candidate_lock_manifest.csv")
standard_model_path <- file.path(project_dir, "data", "processed", "locked_elastic_net_model.rds")
bundle_path <- file.path(project_dir, "data", "processed", "satsa_external_validation_bundle.rds")
checkpoint_path <- file.path(project_dir, "data", "processed", "satsa_beta_extraction_checkpoint.rds")
file_audit_path <- file.path(project_dir, "outputs", "satsa_beta_file_audit.csv")
feature_audit_path <- file.path(project_dir, "outputs", "satsa_feature_coverage.csv")
summary_path <- file.path(project_dir, "outputs", "satsa_beta_extraction_summary.csv")

required_manifest_columns <- c(
  "dataset_id", "resource_id", "resource_type", "remote_url", "local_path",
  "expected_size_mb", "required_for_stage", "description"
)

unquote_fields <- function(values) {
  sub('"$', "", sub('^"', "", values))
}

write_audits <- function(file_audit, feature_audit, summary) {
  dir.create(dirname(file_audit_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(file_audit, file_audit_path, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(feature_audit, feature_audit_path, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(summary, summary_path, row.names = FALSE, fileEncoding = "UTF-8")
}

download_one <- function(url, destination, expected_bytes) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(destination)) {
    actual <- file.info(destination)$size
    if (identical(as.numeric(actual), as.numeric(expected_bytes))) {
      message("이미 완료된 원본 재사용: ", basename(destination))
      return(invisible(destination))
    }
    stop("완성 경로의 파일 크기가 공식 값과 다릅니다: ", destination,
         "\n실제 ", actual, " / 예상 ", expected_bytes,
         "\n파일을 직접 점검한 뒤 잘못된 파일만 제거하세요.")
  }

  partial <- paste0(destination, ".part")
  if (file.exists(partial) && file.info(partial)$size > expected_bytes) {
    stop("부분 다운로드 파일이 공식 크기보다 큽니다: ", partial)
  }
  curl <- Sys.which("curl")
  if (!nzchar(curl)) stop("재개 다운로드에 필요한 curl 실행 파일을 찾지 못했습니다.")
  message("다운로드 시작 또는 재개: ", basename(destination))
  status <- system2(
    curl,
    c(
      "--fail", "--location", "--continue-at", "-", "--retry", "8",
      "--retry-delay", "5", "--output", shQuote(partial), shQuote(url)
    )
  )
  if (!identical(status, 0L)) stop("다운로드가 실패했습니다: ", basename(destination))
  actual <- file.info(partial)$size
  if (!identical(as.numeric(actual), as.numeric(expected_bytes))) {
    stop("다운로드 크기가 공식 값과 다릅니다: ", basename(destination),
         " / 실제 ", actual, " / 예상 ", expected_bytes)
  }
  if (!file.rename(partial, destination)) stop("완료 파일 이름 변경에 실패했습니다: ", destination)
  invisible(destination)
}

read_source_file <- function(path, selected_ids, target_features, extracted, found) {
  connection <- file(path, open = "rt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  header_line <- readLines(connection, n = 1L, warn = FALSE)
  if (length(header_line) != 1L) stop("원본 header를 읽지 못했습니다: ", path)
  header <- unquote_fields(strsplit(header_line, "\t", fixed = TRUE)[[1]])
  if (length(header) != expected_source_samples || anyDuplicated(header)) {
    stop("원본 SID header의 열 수 또는 고유성이 예상과 다릅니다: ", basename(path))
  }
  selected_positions <- match(selected_ids, header)
  if (anyNA(selected_positions)) {
    stop("잠근 SATSA SID가 원본 header에 없습니다: ",
         paste(head(selected_ids[is.na(selected_positions)], 10L), collapse = ", "))
  }

  source_rows <- 0L
  target_rows <- 0L
  repeat {
    lines <- readLines(connection, n = chunk_lines, warn = FALSE)
    if (!length(lines)) break
    source_rows <- source_rows + length(lines)
    first_tabs <- regexpr("\t", lines, fixed = TRUE)
    if (any(first_tabs < 2L)) stop("탭 구분 형식이 아닌 행이 있습니다: ", basename(path))
    probe_ids <- unquote_fields(substring(lines, 1L, first_tabs - 1L))
    target_indices <- match(probe_ids, target_features, nomatch = 0L)
    matched_lines <- which(target_indices > 0L)
    if (!length(matched_lines)) next

    for (line_index in matched_lines) {
      feature_index <- target_indices[[line_index]]
      if (found[[feature_index]]) stop("모델 CpG가 원본 파일 사이에서 중복됩니다: ", target_features[[feature_index]])
      fields <- strsplit(lines[[line_index]], "\t", fixed = TRUE)[[1]]
      if (length(fields) != expected_source_samples + 1L) {
        stop("beta 행의 열 수가 예상과 다릅니다: ", target_features[[feature_index]])
      }
      values <- suppressWarnings(as.numeric(unquote_fields(fields[selected_positions + 1L])))
      nonmissing_text <- nzchar(unquote_fields(fields[selected_positions + 1L])) &
        !tolower(unquote_fields(fields[selected_positions + 1L])) %in% c("na", "nan", "null")
      if (any(is.na(values) & nonmissing_text)) {
        stop("숫자로 변환할 수 없는 beta 값이 있습니다: ", target_features[[feature_index]])
      }
      if (any(values < 0 | values > 1, na.rm = TRUE)) {
        stop("0~1 범위를 벗어난 beta 값이 있습니다: ", target_features[[feature_index]])
      }
      extracted[, feature_index] <- values
      found[[feature_index]] <- TRUE
      target_rows <- target_rows + 1L
    }
  }
  list(
    extracted = extracted, found = found, source_rows = source_rows,
    target_rows = target_rows, header = header
  )
}

for (path in c(download_manifest_path, sample_manifest_path, candidate_model_path,
               candidate_lock_path, standard_model_path)) {
  if (!file.exists(path)) stop("필수 입력 파일이 없습니다: ", path)
}

download_manifest <- read.csv(download_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
missing_columns <- setdiff(required_manifest_columns, names(download_manifest))
if (length(missing_columns)) stop("다운로드 manifest 필수 열이 없습니다: ", paste(missing_columns, collapse = ", "))
resources <- download_manifest[
  download_manifest$dataset_id == accession & download_manifest$resource_type == "beta_matrix",
  , drop = FALSE
]
resources <- resources[order(resources$resource_id), , drop = FALSE]
if (nrow(resources) != 5L || anyDuplicated(resources$resource_id)) {
  stop("SATSA beta 원본은 고유한 5개 resource여야 합니다.")
}
resources$expected_bytes <- round(as.numeric(resources$expected_size_mb) * 1e6)
resources$absolute_path <- file.path(project_dir, resources$local_path)
if (any(!is.finite(resources$expected_bytes)) || any(resources$expected_bytes <= 0)) {
  stop("SATSA 공식 파일 크기가 올바르지 않습니다.")
}

sample_manifest <- read.csv(sample_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_sample_columns <- c(
  "external_sample_id", "participant_id", "family_cluster_id", "chronological_age",
  "sex", "platform", "selection_rule"
)
missing_sample_columns <- setdiff(required_sample_columns, names(sample_manifest))
if (length(missing_sample_columns)) stop("표본 manifest 필수 열이 없습니다: ", paste(missing_sample_columns, collapse = ", "))
if (nrow(sample_manifest) != expected_selected_samples ||
    anyDuplicated(sample_manifest$external_sample_id) || anyDuplicated(sample_manifest$participant_id)) {
  stop("SATSA 잠금 표본은 SID·IID가 고유한 447명이어야 합니다.")
}
if (any(tolower(sample_manifest$platform) != "450k") ||
    any(sample_manifest$selection_rule != "earliest_450k_sample_per_IID")) {
  stop("SATSA 플랫폼 또는 표본 선택 규칙이 잠금 정책과 다릅니다.")
}

candidate_model <- readRDS(candidate_model_path)
standard_model <- readRDS(standard_model_path)
target_features <- as.character(candidate_model$feature_names)
if (length(target_features) != expected_model_features || anyDuplicated(target_features)) {
  stop("후보 모델 feature는 고유한 869개 CpG여야 합니다.")
}
if (!identical(target_features, as.character(standard_model$feature_names))) {
  stop("후보 모델과 표준 모델의 CpG 이름 또는 순서가 다릅니다.")
}
candidate_lock <- read.csv(candidate_lock_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(candidate_lock) != 1L ||
    unname(tools::md5sum(candidate_model_path)) != candidate_lock$model_md5[[1]]) {
  stop("후보 모델과 잠금 manifest의 MD5가 일치하지 않습니다.")
}

present <- file.exists(resources$absolute_path)
actual_bytes <- rep(NA_real_, nrow(resources))
actual_bytes[present] <- file.info(resources$absolute_path[present])$size
size_ok <- present & actual_bytes == resources$expected_bytes
partial_paths <- paste0(resources$absolute_path, ".part")
partial_present <- file.exists(partial_paths)
partial_bytes <- rep(0, nrow(resources))
partial_bytes[partial_present] <- file.info(partial_paths[partial_present])$size
downloaded_bytes <- ifelse(size_ok, resources$expected_bytes, partial_bytes)

cat("외부 후보:", accession, "\n")
cat("공식 원본:", nrow(resources), "개 /", round(sum(resources$expected_bytes) / 1e9, 3), "GB\n")
cat("잠금 표본:", nrow(sample_manifest), "명 / 모델 CpG:", length(target_features), "개\n")
cat("품질 게이트: CpG coverage >=", minimum_feature_coverage,
    "/ 표본별 결측률 <=", maximum_sample_missing_rate, "\n")
print(data.frame(
  resource_id = resources$resource_id,
  expected_bytes = resources$expected_bytes,
  downloaded = present,
  size_ok = size_ok,
  partial_bytes = partial_bytes,
  progress_percent = round(100 * downloaded_bytes / resources$expected_bytes, 2),
  stringsAsFactors = FALSE
), row.names = FALSE)

if (dry_run) {
  cat("dry-run: 원본을 내려받거나 beta를 읽지 않았습니다.\n")
  quit(status = 0)
}

if (download) {
  for (index in seq_len(nrow(resources))) {
    download_one(
      resources$remote_url[[index]], resources$absolute_path[[index]],
      resources$expected_bytes[[index]]
    )
  }
}

if (!extract) {
  cat("추출을 수행하지 않았습니다. --extract --confirm-extraction을 함께 사용하세요.\n")
  quit(status = 0)
}
if (!confirm_extraction) stop("beta 추출에는 --confirm-extraction을 명시해야 합니다. 성능 평가는 수행하지 않습니다.")
if (!all(file.exists(resources$absolute_path))) stop("SATSA beta 5개 파일이 모두 필요합니다. 먼저 --download를 실행하세요.")
actual_bytes <- file.info(resources$absolute_path)$size
if (any(actual_bytes != resources$expected_bytes)) stop("원본 파일 크기가 공식 manifest와 다릅니다.")
if (!force && file.exists(bundle_path)) {
  stop("SATSA 외부 검증 bundle이 이미 있습니다: ", bundle_path,
       "\n재생성 근거가 있을 때만 --force를 사용하세요.")
}

message("원본 5개 파일 MD5 계산 중...")
source_md5 <- unname(tools::md5sum(resources$absolute_path))
signature <- list(
  candidate_model_md5 = unname(tools::md5sum(candidate_model_path)),
  standard_model_md5 = unname(tools::md5sum(standard_model_path)),
  sample_manifest_md5 = unname(tools::md5sum(sample_manifest_path)),
  resource_ids = resources$resource_id,
  source_md5 = source_md5,
  expected_bytes = resources$expected_bytes,
  target_features = target_features,
  selected_ids = as.character(sample_manifest$external_sample_id)
)

if (force && file.exists(checkpoint_path)) unlink(checkpoint_path)
checkpoint <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else list(
  signature = signature,
  extracted = matrix(
    NA_real_, nrow = nrow(sample_manifest), ncol = length(target_features),
    dimnames = list(sample_manifest$external_sample_id, target_features)
  ),
  found = rep(FALSE, length(target_features)),
  completed = character(),
  file_audit = list(),
  total_source_rows = 0L
)
if (!identical(checkpoint$signature, signature)) {
  stop("기존 추출 checkpoint의 입력 해시 또는 잠금 목록이 현재 실행과 다릅니다.")
}

for (index in seq_len(nrow(resources))) {
  resource_id <- resources$resource_id[[index]]
  if (resource_id %in% checkpoint$completed) {
    message("checkpoint 재사용: ", resource_id)
    next
  }
  message("스트리밍 추출 시작: ", resource_id)
  result <- read_source_file(
    resources$absolute_path[[index]], signature$selected_ids, target_features,
    checkpoint$extracted, checkpoint$found
  )
  checkpoint$extracted <- result$extracted
  checkpoint$found <- result$found
  checkpoint$total_source_rows <- checkpoint$total_source_rows + result$source_rows
  checkpoint$file_audit[[resource_id]] <- data.frame(
    accession = accession,
    resource_id = resource_id,
    local_path = resources$local_path[[index]],
    expected_bytes = resources$expected_bytes[[index]],
    actual_bytes = actual_bytes[[index]],
    size_ok = actual_bytes[[index]] == resources$expected_bytes[[index]],
    md5 = source_md5[[index]],
    source_rows = result$source_rows,
    source_samples = length(result$header),
    selected_sid_matches = sum(signature$selected_ids %in% result$header),
    extracted_target_rows = result$target_rows,
    stringsAsFactors = FALSE
  )
  checkpoint$completed <- c(checkpoint$completed, resource_id)
  saveRDS(checkpoint, checkpoint_path)
  message("완료: ", resource_id, " / source rows=", result$source_rows,
          " / 모델 CpG=", result$target_rows)
}

file_audit <- do.call(rbind, checkpoint$file_audit[resources$resource_id])
row.names(file_audit) <- NULL
if (checkpoint$total_source_rows != expected_source_rows) {
  stop("SATSA 원본 CpG 총 행 수가 공식 README의 250,816개와 다릅니다: ", checkpoint$total_source_rows)
}

feature_audit <- data.frame(
  feature = target_features,
  available = checkpoint$found,
  missing_all_samples = colSums(is.na(checkpoint$extracted)) == nrow(checkpoint$extracted),
  selected_sample_missing_rate = colMeans(is.na(checkpoint$extracted)),
  stringsAsFactors = FALSE
)
feature_coverage <- mean(feature_audit$available)
sample_missing_rate <- rowMeans(is.na(checkpoint$extracted))
values <- checkpoint$extracted[is.finite(checkpoint$extracted)]
range_ok <- length(values) > 0L && min(values) >= 0 && max(values) <= 1
quality_pass <- feature_coverage >= minimum_feature_coverage &&
  max(sample_missing_rate) <= maximum_sample_missing_rate && range_ok

summary <- data.frame(
  accession = accession,
  source_files = nrow(resources),
  source_bytes = sum(actual_bytes),
  source_rows = checkpoint$total_source_rows,
  source_samples = expected_source_samples,
  selected_samples = nrow(sample_manifest),
  model_features = length(target_features),
  available_features = sum(feature_audit$available),
  missing_features = sum(!feature_audit$available),
  feature_coverage = feature_coverage,
  sample_missing_rate_max = max(sample_missing_rate),
  sample_missing_rate_mean = mean(sample_missing_rate),
  beta_min = if (length(values)) min(values) else NA_real_,
  beta_max = if (length(values)) max(values) else NA_real_,
  duplicate_selected_sid = anyDuplicated(signature$selected_ids),
  total_sid_header_match = all(file_audit$selected_sid_matches == expected_selected_samples),
  quality_gate = if (quality_pass) "PASS" else "FAIL",
  external_evaluation_used = FALSE,
  candidate_model_md5 = signature$candidate_model_md5,
  standard_model_md5 = signature$standard_model_md5,
  sample_manifest_md5 = signature$sample_manifest_md5,
  extracted_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
write_audits(file_audit, feature_audit, summary)

message("=== SATSA beta 추출 품질검사 ===")
print(summary, row.names = FALSE)
if (!quality_pass) {
  stop("SATSA beta가 사전 품질 게이트를 통과하지 못했습니다. bundle과 성능 결과를 만들지 않았습니다.")
}

bundle <- list(
  x = checkpoint$extracted,
  y = as.numeric(sample_manifest$chronological_age),
  sample_id = as.character(sample_manifest$external_sample_id),
  participant_id = as.character(sample_manifest$participant_id),
  family_cluster_id = as.character(sample_manifest$family_cluster_id),
  sex = as.character(sample_manifest$sex),
  feature_names = target_features,
  source_accession = accession,
  source_file_md5 = stats::setNames(source_md5, resources$resource_id),
  sample_manifest_md5 = signature$sample_manifest_md5,
  candidate_model_md5 = signature$candidate_model_md5,
  standard_model_md5 = signature$standard_model_md5,
  selection_rule = "earliest_450k_sample_per_IID",
  external_evaluation_used = FALSE
)
dir.create(dirname(bundle_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(bundle, bundle_path)
message("외부 검증 bundle 저장 완료: ", bundle_path)
message("성능 평가는 수행하지 않았습니다.")
