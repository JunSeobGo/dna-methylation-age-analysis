#!/usr/bin/env Rscript

# 1차 학습 코호트(model_policy == include) 8개를 모아 Elastic Net 중첩 교차검증에
# 바로 투입할 수 있는 학습용 bundle을 생성한다.
#
# 원칙
# - 869개 공통 CpG만 동일한 순서로 추출한다.
# - 표본 방향(행=표본, 열=CpG)으로 정렬한다.
# - 결측치(NA)는 대치하지 않고 bundle에 그대로 보존한다.
#   결측 대치와 표준화는 교차검증의 각 학습 fold 안에서만 수행해야 하므로
#   이 단계에서 전체 데이터 기준 대치를 하지 않는다.
# - 원본과 산출물은 Git에 포함하지 않는다.

args <- commandArgs(trailingOnly = TRUE)
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

manifest_path <- file.path(project_dir, "config", "gse207605_cohort_manifest.csv")
raw_dir <- file.path(project_dir, "data", "raw", "GSE207605")
common_probe_path <- file.path(project_dir, "data", "processed", "gse207605_primary_common_probes.txt")
bundle_path <- file.path(project_dir, "data", "processed", "gse207605_training_bundle.rds")
summary_path <- file.path(project_dir, "outputs", "gse207605_training_bundle_summary.csv")

# --- 입력 점검 ---------------------------------------------------------------

if (!file.exists(manifest_path)) stop("코호트 명세를 찾을 수 없습니다: ", manifest_path)
if (!file.exists(common_probe_path)) stop("공통 CpG 목록을 찾을 수 없습니다: ", common_probe_path)

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("source_series_id", "expected_samples", "model_policy")
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns)) {
  stop("코호트 명세에 필요한 열이 없습니다: ", paste(missing_columns, collapse = ", "))
}
if (!all(manifest$model_policy %in% c("include", "conditional", "exclude"))) {
  stop("model_policy는 include, conditional, exclude 중 하나여야 합니다.")
}

include_manifest <- manifest[manifest$model_policy == "include", , drop = FALSE]
include_manifest <- include_manifest[order(include_manifest$source_series_id), , drop = FALSE]
if (nrow(include_manifest) != 8L) {
  stop("1차 학습 코호트는 8개여야 하지만 ", nrow(include_manifest), "개가 확인되었습니다.")
}

common_probes <- readLines(common_probe_path, warn = FALSE)
common_probes <- common_probes[nzchar(common_probes)]
if (anyDuplicated(common_probes)) stop("공통 CpG 목록에 중복이 있습니다.")
if (length(common_probes) != 869L) {
  stop("공통 CpG는 869개여야 하지만 ", length(common_probes), "개가 확인되었습니다.")
}

expected_total <- sum(include_manifest$expected_samples)

message("1차 학습 코호트 ", nrow(include_manifest), "개, 예상 표본 ", expected_total, "명")
message("공통 CpG ", length(common_probes), "개를 동일 순서로 추출합니다.")

# --- 코호트별 추출 -----------------------------------------------------------

extract_one_cohort <- function(resource) {
  series_id <- resource$source_series_id[[1]]
  path <- file.path(raw_dir, paste0("GSE207605_", series_id, ".csv.gz"))
  if (!file.exists(path)) stop("원본 파일이 없습니다: ", path)

  data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(data) < 2L || nrow(data) < 2L) stop("행렬 구조가 올바르지 않습니다: ", series_id)
  if (as.character(data[[1]][[1]]) != "age") stop("첫 데이터 행이 age가 아닙니다: ", series_id)

  sample_ids <- names(data)[-1]
  ages <- suppressWarnings(as.numeric(data[1, -1, drop = TRUE]))
  probes <- as.character(data[-1, 1, drop = TRUE])
  beta <- as.matrix(data[-1, -1, drop = FALSE])
  storage.mode(beta) <- "numeric"

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

  # 869개 공통 CpG를 목록 순서 그대로 찾는다. 하나라도 없으면 공통 집합 정의와 모순이므로 중단한다.
  probe_index <- match(common_probes, probes)
  if (any(is.na(probe_index))) {
    missing_probes <- common_probes[is.na(probe_index)]
    stop(
      series_id, " 코호트에서 공통 CpG ", length(missing_probes),
      "개를 찾을 수 없습니다 (예: ", paste(head(missing_probes, 3), collapse = ", "), ")"
    )
  }

  # beta는 (CpG × 표본)이므로 공통 CpG만 고른 뒤 전치하여 (표본 × CpG)로 만든다.
  beta_common <- beta[probe_index, , drop = FALSE]
  x <- t(beta_common)
  rownames(x) <- sample_ids
  colnames(x) <- common_probes

  list(
    series_id = series_id,
    x = x,
    age = ages,
    sample_ids = sample_ids
  )
}

results <- lapply(seq_len(nrow(include_manifest)), function(index) {
  extract_one_cohort(include_manifest[index, , drop = FALSE])
})

# --- 통합 및 정렬 ------------------------------------------------------------

x <- do.call(rbind, lapply(results, `[[`, "x"))
y <- unlist(lapply(results, `[[`, "age"), use.names = FALSE)
sample_id <- unlist(lapply(results, `[[`, "sample_ids"), use.names = FALSE)
group <- unlist(
  lapply(results, function(item) rep(item$series_id, length(item$sample_ids))),
  use.names = FALSE
)
feature_names <- common_probes

# --- 누수 방지 및 무결성 검증 ------------------------------------------------

stopifnot(
  "X 행 수가 예상 표본 수와 다릅니다" = nrow(x) == expected_total,
  "X 열 수가 869가 아닙니다" = ncol(x) == 869L,
  "y 길이가 표본 수와 다릅니다" = length(y) == expected_total,
  "group 길이가 표본 수와 다릅니다" = length(group) == expected_total,
  "sample_id 길이가 표본 수와 다릅니다" = length(sample_id) == expected_total,
  "feature_names 길이가 869가 아닙니다" = length(feature_names) == 869L
)

# 행 순서와 메타데이터 순서가 정확히 일치하는지 확인한다.
if (!identical(rownames(x), sample_id)) stop("X 행 이름과 sample_id 순서가 일치하지 않습니다.")
if (!identical(colnames(x), feature_names)) stop("X 열 이름과 feature_names 순서가 일치하지 않습니다.")

if (anyDuplicated(sample_id)) stop("코호트 사이에 중복 GSM ID가 있습니다.")
if (anyDuplicated(feature_names)) stop("feature_names에 중복 CpG가 있습니다.")

# 코호트별 표본 수가 명세와 일치하는지 확인한다.
observed_counts <- table(group)[include_manifest$source_series_id]
expected_counts <- setNames(include_manifest$expected_samples, include_manifest$source_series_id)
if (!all(observed_counts == expected_counts)) {
  stop("코호트별 표본 수가 명세와 일치하지 않습니다.")
}

# beta 범위 검증: NA를 제외한 값이 모두 0~1 안에 있어야 한다.
finite_x <- x[is.finite(x)]
if (length(finite_x) && (min(finite_x) < 0 || max(finite_x) > 1)) {
  stop("0~1 범위를 벗어난 beta 값이 있습니다: min=", min(finite_x), ", max=", max(finite_x))
}

# 결측치 보고 (대치하지 않는다).
na_total <- sum(is.na(x))
na_ratio <- na_total / length(x)
na_by_cohort <- sapply(split(seq_len(nrow(x)), group), function(rows) sum(is.na(x[rows, , drop = FALSE])))
na_by_cohort <- na_by_cohort[include_manifest$source_series_id]

message("결측치 총 ", na_total, "개 (비율 ", signif(na_ratio, 4), ")")

# --- bundle 저장 -------------------------------------------------------------

bundle <- list(
  x = x,
  y = y,
  group = group,
  sample_id = sample_id,
  feature_names = feature_names
)

dir.create(dirname(bundle_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(bundle, bundle_path)

# --- 생성 요약(outputs/) -----------------------------------------------------

per_cohort <- do.call(rbind, lapply(include_manifest$source_series_id, function(series_id) {
  rows <- which(group == series_id)
  ages <- y[rows]
  block <- x[rows, , drop = FALSE]
  data.frame(
    source_series_id = series_id,
    samples = length(rows),
    expected_samples = expected_counts[[series_id]],
    age_min = min(ages),
    age_median = median(ages),
    age_max = max(ages),
    na_cells = na_by_cohort[[series_id]],
    na_ratio = na_by_cohort[[series_id]] / length(block),
    stringsAsFactors = FALSE
  )
}))

overall <- data.frame(
  source_series_id = "ALL",
  samples = nrow(x),
  expected_samples = expected_total,
  age_min = min(y),
  age_median = median(y),
  age_max = max(y),
  na_cells = na_total,
  na_ratio = na_ratio,
  stringsAsFactors = FALSE
)

summary_table <- rbind(overall, per_cohort)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, summary_path, row.names = FALSE, fileEncoding = "UTF-8")

# --- 최종 보고 ---------------------------------------------------------------

message("학습 bundle 저장 완료: ", bundle_path)
message(sprintf("X: %d행 × %d열, 코호트 %d개", nrow(x), ncol(x), length(unique(group))))
print(summary_table, row.names = FALSE)
message("생성 요약 저장: ", summary_path)
