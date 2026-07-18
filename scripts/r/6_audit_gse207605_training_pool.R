#!/usr/bin/env Rscript

# GSE207605의 19개 코호트를 내려받고 표본 수, 연령, 공통 CpG, 코호트 간 중복을 감사한다.
# 질환 혼합 및 반복측정 가능성이 있는 코호트는 명세의 정책에 따라 1차 학습 풀에서 분리한다.

args <- commandArgs(trailingOnly = TRUE)
download_missing <- "--download" %in% args
dry_run <- "--dry-run" %in% args
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

manifest_path <- file.path(project_dir, "config", "gse207605_cohort_manifest.csv")
raw_dir <- file.path(project_dir, "data", "raw", "GSE207605")
audit_path <- file.path(project_dir, "outputs", "gse207605_cohort_audit.csv")
summary_path <- file.path(project_dir, "outputs", "gse207605_policy_summary.csv")
sample_manifest_path <- file.path(project_dir, "data", "processed", "gse207605_sample_manifest.csv")
common_probe_path <- file.path(project_dir, "data", "processed", "gse207605_primary_common_probes.txt")

if (!file.exists(manifest_path)) stop("코호트 명세를 찾을 수 없습니다: ", manifest_path)

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c(
  "source_series_id", "expected_samples", "proposed_role", "model_policy",
  "reason", "remote_url", "source_url"
)
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns)) {
  stop("코호트 명세에 필요한 열이 없습니다: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(manifest$source_series_id)) stop("코호트 명세에 중복 series ID가 있습니다.")
if (!all(manifest$model_policy %in% c("include", "conditional", "exclude"))) {
  stop("model_policy는 include, conditional, exclude 중 하나여야 합니다.")
}
if (nrow(manifest) != 19L) stop("GSE207605 구성 코호트는 19개여야 합니다.")

message("GSE207605 코호트 정책")
print(manifest[, c("source_series_id", "expected_samples", "model_policy", "proposed_role")], row.names = FALSE)
if (dry_run) quit(status = 0)

is_gzip_file <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  identical(readBin(connection, what = "raw", n = 2L), as.raw(c(0x1f, 0x8b)))
}

download_resource <- function(series_id, remote_url) {
  destination <- file.path(raw_dir, paste0("GSE207605_", series_id, ".csv.gz"))
  partial_destination <- paste0(destination, ".part")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(destination)) {
    if (!is_gzip_file(destination)) stop("gzip 형식이 아닌 기존 파일입니다: ", destination)
    message("기존 파일 사용: ", series_id)
    return(destination)
  }
  if (!download_missing) {
    stop("원본 파일이 없습니다: ", destination, "\n--download 옵션으로 내려받으세요.")
  }

  curl_path <- Sys.which("curl.exe")
  if (!nzchar(curl_path)) stop("재개 다운로드에 필요한 curl.exe를 찾을 수 없습니다.")
  message("다운로드: ", series_id)
  status <- system2(
    curl_path,
    c(
      "--fail", "--location", "--continue-at", "-", "--retry", "5",
      "--retry-delay", "3", "--output", partial_destination, remote_url
    )
  )
  if (!identical(status, 0L)) stop("다운로드에 실패했습니다: ", series_id)
  if (!is_gzip_file(partial_destination)) stop("다운로드 파일의 gzip 검증에 실패했습니다: ", series_id)
  if (!file.rename(partial_destination, destination)) stop("완료 파일 이름 변경에 실패했습니다: ", series_id)
  destination
}

audit_one_cohort <- function(resource) {
  series_id <- resource$source_series_id[[1]]
  path <- download_resource(series_id, resource$remote_url[[1]])
  data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  if (ncol(data) < 2L || nrow(data) < 2L) stop("행렬 구조가 올바르지 않습니다: ", series_id)
  if (as.character(data[[1]][[1]]) != "age") stop("첫 데이터 행이 age가 아닙니다: ", series_id)

  sample_ids <- names(data)[-1]
  ages <- suppressWarnings(as.numeric(data[1, -1, drop = TRUE]))
  probes <- as.character(data[-1, 1, drop = TRUE])
  beta <- as.matrix(data[-1, -1, drop = FALSE])
  storage.mode(beta) <- "numeric"
  finite_beta <- is.finite(beta)

  observed_samples <- length(sample_ids)
  if (observed_samples != resource$expected_samples[[1]]) {
    stop(
      series_id, " 표본 수가 명세와 다릅니다: expected=", resource$expected_samples[[1]],
      ", observed=", observed_samples
    )
  }
  if (anyDuplicated(sample_ids)) stop("코호트 내부에 중복 GSM ID가 있습니다: ", series_id)
  if (anyDuplicated(probes)) stop("코호트 내부에 중복 CpG ID가 있습니다: ", series_id)
  if (any(!is.finite(ages))) stop("연령 결측 또는 비수치 값이 있습니다: ", series_id)

  sample_manifest <- data.frame(
    source_series_id = series_id,
    sample_id = sample_ids,
    chronological_age = ages,
    model_policy = resource$model_policy[[1]],
    proposed_role = resource$proposed_role[[1]],
    stringsAsFactors = FALSE
  )
  audit <- data.frame(
    source_series_id = series_id,
    model_policy = resource$model_policy[[1]],
    samples = observed_samples,
    unique_ages = length(unique(ages)),
    age_min = min(ages),
    age_median = median(ages),
    age_max = max(ages),
    probe_count = length(probes),
    beta_missing = sum(!finite_beta),
    beta_out_of_range = sum(beta[finite_beta] < 0 | beta[finite_beta] > 1),
    stringsAsFactors = FALSE
  )
  list(audit = audit, sample_manifest = sample_manifest, probes = probes)
}

results <- lapply(seq_len(nrow(manifest)), function(index) {
  audit_one_cohort(manifest[index, , drop = FALSE])
})

audit <- do.call(rbind, lapply(results, `[[`, "audit"))
sample_manifest <- do.call(rbind, lapply(results, `[[`, "sample_manifest"))
all_sample_ids <- sample_manifest$sample_id
if (anyDuplicated(all_sample_ids)) {
  duplicated_ids <- unique(all_sample_ids[duplicated(all_sample_ids)])
  stop("코호트 사이에 중복 GSM ID가 있습니다: ", paste(duplicated_ids, collapse = ", "))
}
if (any(audit$beta_out_of_range > 0L)) stop("0~1 범위를 벗어난 beta 값이 있습니다.")

include_ids <- manifest$source_series_id[manifest$model_policy == "include"]
probe_sets <- setNames(lapply(results, `[[`, "probes"), manifest$source_series_id)
common_probes <- sort(Reduce(intersect, probe_sets[include_ids]))
if (!length(common_probes)) stop("1차 학습 코호트 사이의 공통 CpG가 없습니다.")

policy_levels <- c("include", "conditional", "exclude")
summary <- do.call(rbind, lapply(policy_levels, function(policy) {
  selected <- audit$model_policy == policy
  data.frame(
    model_policy = policy,
    cohorts = sum(selected),
    samples = sum(audit$samples[selected]),
    stringsAsFactors = FALSE
  )
}))
summary$primary_common_probes <- ifelse(summary$model_policy == "include", length(common_probes), NA_integer_)

dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(sample_manifest_path), recursive = TRUE, showWarnings = FALSE)
write.csv(audit, audit_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(summary, summary_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sample_manifest, sample_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")
writeLines(common_probes, common_probe_path, useBytes = TRUE)

message("감사 완료: 총 ", nrow(sample_manifest), "개 표본, 코호트 간 중복 GSM 0개")
print(audit, row.names = FALSE)
print(summary, row.names = FALSE)
message("1차 학습 공통 CpG: ", length(common_probes), "개")
