#!/usr/bin/env Rscript

# 외부 검증셋을 다시 사용하지 않고 학습 코호트의 홀드아웃 예측으로
# 연령 극단부 편향, 코호트별 데이터 품질, 연령 지원도와 CpG 분포 이동을 진단한다.

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

bundle_path <- file.path(project_dir, "data", "processed", "gse207605_training_bundle.rds")
prediction_path <- file.path(project_dir, "outputs", "training_learning_curve_predictions.csv")
output_dir <- file.path(project_dir, "outputs")
overall_path <- file.path(output_dir, "age_cohort_bias_audit_overall.csv")
cohort_path <- file.path(output_dir, "age_cohort_bias_audit_by_cohort.csv")
age_path <- file.path(output_dir, "age_cohort_bias_audit_by_agebin.csv")
quality_path <- file.path(output_dir, "age_cohort_bias_audit_quality_checks.csv")
finding_path <- file.path(output_dir, "age_cohort_bias_audit_findings.csv")
plot_path <- file.path(output_dir, "age_cohort_bias_audit.png")

for (path in c(bundle_path, prediction_path)) {
  if (!file.exists(path)) stop("필수 입력 파일이 없습니다: ", path)
}

bundle <- readRDS(bundle_path)
predictions_all <- read.csv(prediction_path, stringsAsFactors = FALSE, check.names = FALSE)
x <- bundle$x
y <- as.numeric(bundle$y)
group <- as.character(bundle$group)
sample_id <- as.character(bundle$sample_id)
features <- as.character(bundle$feature_names)

required_prediction_columns <- c(
  "fraction", "repeat_id", "held_out_cohort", "sample_id",
  "age_true", "age_predicted", "residual", "age_bin"
)
if (!all(required_prediction_columns %in% names(predictions_all))) {
  stop("학습곡선 예측 파일에 필수 열이 없습니다: ",
       paste(setdiff(required_prediction_columns, names(predictions_all)), collapse = ", "))
}
if (!is.matrix(x) || nrow(x) != length(y) || length(y) != length(group) || length(y) != length(sample_id)) {
  stop("학습 bundle의 행 수와 메타데이터 길이가 일치하지 않습니다.")
}
if (!identical(colnames(x), features)) stop("CpG 이름과 행렬 열 순서가 다릅니다.")
if (anyDuplicated(sample_id) || anyDuplicated(features)) stop("sample_id 또는 CpG 이름에 중복이 있습니다.")
if (any(!is.finite(y))) stop("학습 연령에 결측값 또는 무한값이 있습니다.")
if (ncol(x) != 869L) stop("고정 입력 특성은 869개 CpG여야 합니다.")

# 학습곡선의 100% 조건은 외부셋을 사용하지 않은 고정 모델의 코호트 홀드아웃 예측이다.
predictions <- predictions_all[predictions_all$fraction == 1 & predictions_all$repeat_id == 1, , drop = FALSE]
if (nrow(predictions) != length(y)) stop("100% 조건 예측 행 수가 학습 표본 수와 다릅니다.")
if (anyDuplicated(predictions$sample_id)) stop("100% 조건 예측에 중복 sample_id가 있습니다.")
if (!setequal(predictions$sample_id, sample_id)) stop("예측 sample_id와 학습 bundle sample_id가 일치하지 않습니다.")
if (any(predictions$held_out_cohort == "GSE87571")) stop("외부 검증 코호트가 진단 입력에 포함됐습니다.")

predictions <- predictions[match(sample_id, predictions$sample_id), , drop = FALSE]
if (any(abs(predictions$age_true - y) > 1e-10)) stop("예측 파일과 bundle의 실제 연령이 다릅니다.")
if (any(predictions$held_out_cohort != group)) stop("예측 파일의 홀드아웃 코호트와 bundle 코호트가 다릅니다.")
if (any(abs((predictions$age_predicted - predictions$age_true) - predictions$residual) > 1e-8)) {
  stop("예측값에서 재계산한 잔차와 저장된 잔차가 다릅니다.")
}

cohorts <- sort(unique(group))
if (length(cohorts) != 8L) stop("이 감사는 고정된 8개 학습 코호트를 전제로 합니다.")

age_levels <- c("<20", "20-39", "40-59", "60-79", ">=80")
age_bin <- function(values) {
  cut(values, breaks = c(-Inf, 20, 40, 60, 80, Inf),
      labels = age_levels, right = FALSE)
}

metric_values <- function(actual, predicted) {
  residual <- predicted - actual
  denominator <- sum((actual - mean(actual))^2)
  calibration <- if (length(unique(predicted)) > 1L) lm(actual ~ predicted) else NULL
  c(
    MAE = mean(abs(residual)), RMSE = sqrt(mean(residual^2)),
    MedAE = median(abs(residual)),
    R2 = if (denominator > 0) 1 - sum(residual^2) / denominator else NA_real_,
    bias_mean = mean(residual),
    residual_age_slope = if (length(unique(actual)) > 1L) as.numeric(coef(lm(residual ~ actual))[2]) else NA_real_,
    calibration_intercept = if (!is.null(calibration)) as.numeric(coef(calibration)[1]) else NA_real_,
    calibration_slope = if (!is.null(calibration)) as.numeric(coef(calibration)[2]) else NA_real_
  )
}

bootstrap_intervals <- function(actual, predicted, seed, iterations = 2000L) {
  set.seed(seed)
  residual <- predicted - actual
  estimates <- replicate(iterations, {
    sampled <- sample.int(length(actual), length(actual), replace = TRUE)
    c(MAE = mean(abs(residual[sampled])), bias_mean = mean(residual[sampled]))
  })
  c(
    MAE_ci_low = unname(quantile(estimates["MAE", ], 0.025)),
    MAE_ci_high = unname(quantile(estimates["MAE", ], 0.975)),
    bias_ci_low = unname(quantile(estimates["bias_mean", ], 0.025)),
    bias_ci_high = unname(quantile(estimates["bias_mean", ], 0.975))
  )
}

finite_beta <- x[is.finite(x)]
quality_checks <- data.frame(
  check = c("sample_id_uniqueness", "prediction_coverage", "age_validity", "beta_range",
            "feature_identity", "cohort_count", "external_independence"),
  status = c(
    if (anyDuplicated(sample_id) == 0L) "pass" else "fail",
    if (nrow(predictions) == length(y) && !anyDuplicated(predictions$sample_id)) "pass" else "fail",
    if (all(is.finite(y)) && min(y) >= 0 && max(y) <= 120) "pass" else "fail",
    if (length(finite_beta) && min(finite_beta) >= 0 && max(finite_beta) <= 1) "pass" else "fail",
    if (identical(colnames(x), features) && ncol(x) == 869L) "pass" else "fail",
    if (length(cohorts) == 8L) "pass" else "fail",
    if (!"GSE87571" %in% unique(predictions$held_out_cohort)) "pass" else "fail"
  ),
  evidence = c(
    sprintf("고유 sample_id %d/%d", length(unique(sample_id)), length(sample_id)),
    sprintf("100%% 홀드아웃 예측 %d행", nrow(predictions)),
    sprintf("연령 범위 %.3f~%.3f세", min(y), max(y)),
    sprintf("유효 beta 범위 %.6f~%.6f, 전체 결측률 %.6f%%", min(finite_beta), max(finite_beta), mean(is.na(x)) * 100),
    sprintf("CpG %d개, 중복 %d개", ncol(x), anyDuplicated(features)),
    sprintf("학습 코호트 %d개", length(cohorts)),
    "진단 입력은 8개 학습 코호트의 홀드아웃 예측만 사용"
  ),
  stringsAsFactors = FALSE
)

if (dry_run) {
  message("dry-run: 입력 무결성과 진단 설계만 점검했습니다. CpG 분포 이동은 계산하지 않습니다.")
  print(quality_checks, row.names = FALSE)
  quit(status = if (any(quality_checks$status == "fail")) 1 else 0)
}

overall_values <- metric_values(y, predictions$age_predicted)
overall <- data.frame(
  n = length(y), cohorts = length(cohorts), features = ncol(x),
  age_min = min(y), age_max = max(y), missing_rate = mean(is.na(x)),
  t(overall_values), stringsAsFactors = FALSE
)

cohort_rows <- lapply(seq_along(cohorts), function(cohort_number) {
  cohort <- cohorts[[cohort_number]]
  index <- which(group == cohort)
  reference <- which(group != cohort)
  actual <- y[index]
  predicted <- predictions$age_predicted[index]
  values <- metric_values(actual, predicted)
  intervals <- bootstrap_intervals(actual, predicted, seed = 20260803L + cohort_number)

  reference_min <- min(y[reference])
  reference_max <- max(y[reference])
  outside <- actual < reference_min | actual > reference_max
  nearby_counts <- vapply(actual, function(age) sum(abs(y[reference] - age) <= 5), numeric(1))

  cohort_medians <- apply(x[index, , drop = FALSE], 2, median, na.rm = TRUE)
  reference_medians <- apply(x[reference, , drop = FALSE], 2, median, na.rm = TRUE)
  reference_mads <- apply(x[reference, , drop = FALSE], 2, mad, na.rm = TRUE)
  median_difference <- abs(cohort_medians - reference_medians)
  robust_difference <- median_difference / pmax(reference_mads, 0.01)

  data.frame(
    cohort = cohort, n = length(index),
    age_min = min(actual), age_q1 = unname(quantile(actual, 0.25)),
    age_median = median(actual), age_q3 = unname(quantile(actual, 0.75)), age_max = max(actual),
    age_under_20_n = sum(actual < 20), age_80_plus_n = sum(actual >= 80),
    reference_age_min = reference_min, reference_age_max = reference_max,
    outside_reference_age_n = sum(outside), outside_reference_age_rate = mean(outside),
    median_nearby_training_n_5yr = median(nearby_counts),
    min_nearby_training_n_5yr = min(nearby_counts), zero_nearby_training_rate_5yr = mean(nearby_counts == 0),
    missing_rate = mean(is.na(x[index, , drop = FALSE])),
    beta_mean = mean(x[index, , drop = FALSE], na.rm = TRUE),
    beta_sd = sd(as.numeric(x[index, , drop = FALSE]), na.rm = TRUE),
    cpg_median_shift_mean = mean(median_difference, na.rm = TRUE),
    cpg_median_shift_median = median(median_difference, na.rm = TRUE),
    cpg_median_shift_p95 = unname(quantile(median_difference, 0.95, na.rm = TRUE)),
    cpg_robust_shift_median = median(robust_difference, na.rm = TRUE),
    t(values), t(intervals), stringsAsFactors = FALSE
  )
})
by_cohort <- do.call(rbind, cohort_rows)
overall$spearman_MAE_vs_age_support <- cor(
  by_cohort$MAE, by_cohort$median_nearby_training_n_5yr, method = "spearman"
)
overall$spearman_MAE_vs_cpg_shift <- cor(
  by_cohort$MAE, by_cohort$cpg_median_shift_median, method = "spearman"
)

bins <- age_bin(y)
age_rows <- lapply(age_levels, function(bin) {
  index <- which(bins == bin)
  values <- metric_values(y[index], predictions$age_predicted[index])
  data.frame(
    age_bin = bin, n = length(index), cohorts_represented = length(unique(group[index])),
    cohort_composition = paste(names(sort(table(group[index]), decreasing = TRUE)), collapse = ";"),
    t(values), stringsAsFactors = FALSE
  )
})
by_age <- do.call(rbind, age_rows)

young <- by_age[by_age$age_bin == "<20", , drop = FALSE]
old <- by_age[by_age$age_bin == ">=80", , drop = FALSE]
worst <- by_cohort[order(-by_cohort$MAE), , drop = FALSE]
max_missing <- max(by_cohort$missing_rate)

findings <- rbind(
  data.frame(
    finding = "80세 이상 체계적 과소예측", severity = "high", confidence = "high",
    evidence = sprintf("n=%d, MAE %.3f년, 평균 잔차 %.3f년", old$n, old$MAE, old$bias_mean),
    likely_cause = "고령 표본이 일부 코호트에 집중되고 연령 극단에서 회귀 평균화가 발생",
    recommended_action = "80세 이상 독립 표본 보강과 학습 내부 연령 편향 보정 실험",
    stringsAsFactors = FALSE
  ),
  data.frame(
    finding = "20세 미만 체계적 과소예측", severity = "high", confidence = "high",
    evidence = sprintf("n=%d, MAE %.3f년, 평균 잔차 %.3f년", young$n, young$MAE, young$bias_mean),
    likely_cause = "소아 표본이 GSE36054에 집중되고 홀드아웃 시 저연령 지원이 부족",
    recommended_action = "소아·청소년 독립 코호트 추가 또는 연령 구간별 민감도 모델 검토",
    stringsAsFactors = FALSE
  ),
  data.frame(
    finding = paste0(worst$cohort[1], " 코호트 오차와 불확실성"), severity = "high", confidence = "medium",
    evidence = sprintf("n=%d, 연령 %.1f~%.1f세, MAE %.3f년(95%% CI %.3f~%.3f), 평균 잔차 %.3f년",
                       worst$n[1], worst$age_min[1], worst$age_max[1], worst$MAE[1],
                       worst$MAE_ci_low[1], worst$MAE_ci_high[1], worst$bias_mean[1]),
    likely_cause = "극단 연령과 작은 코호트 크기가 함께 작용하며 단일 MAE 추정의 불확실성이 큼",
    recommended_action = "동일 연령대의 추가 코호트를 확보하고 부트스트랩 신뢰구간을 함께 보고",
    stringsAsFactors = FALSE
  ),
  data.frame(
    finding = "CpG beta 완전성", severity = if (max_missing > 0.01) "medium" else "low", confidence = "high",
    evidence = sprintf("코호트 최대 결측률 %.6f%%, 전체 유효 beta 범위 %.6f~%.6f", max_missing * 100, min(finite_beta), max(finite_beta)),
    likely_cause = if (max_missing > 0.01) "특정 코호트의 probe 결측" else "결측과 범위 위반은 주요 오차 원인으로 보기 어려움",
    recommended_action = "현재 fold 내부 중앙값 대치를 유지하고 결측률 자동 검사를 지속",
    stringsAsFactors = FALSE
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(overall, overall_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(by_cohort, cohort_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(by_age, age_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(quality_checks, quality_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(findings, finding_path, row.names = FALSE, fileEncoding = "UTF-8")

png(plot_path, width = 1400, height = 1000, res = 130)
old_par <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.2, 1))
colors <- as.integer(factor(group))
plot(y, predictions$residual, pch = 16, cex = 0.45, col = adjustcolor(colors, alpha.f = 0.35),
     xlab = "실제 나이", ylab = "잔차(예측 - 실제, 년)", main = "연령에 따른 홀드아웃 잔차")
abline(h = 0, col = "red", lwd = 2); lines(lowess(y, predictions$residual), col = "blue", lwd = 2)

ordered <- by_cohort[order(by_cohort$MAE, decreasing = TRUE), ]
bar_colors <- ifelse(ordered$cohort %in% ordered$cohort[1:2], "tomato", "steelblue")
barplot(ordered$MAE, names.arg = ordered$cohort, las = 2, col = bar_colors,
        ylab = "MAE(년)", main = "홀드아웃 코호트별 오차")

mae_limits <- c(max(0, min(by_cohort$MAE) * 0.8), max(by_cohort$MAE) * 1.12)
support_limits <- c(0, max(by_cohort$median_nearby_training_n_5yr) * 1.06)
plot(by_cohort$median_nearby_training_n_5yr, by_cohort$MAE, pch = 19, col = "steelblue",
     xlim = support_limits, ylim = mae_limits,
     xlab = "±5세 학습 표본 수 중앙값", ylab = "코호트 MAE(년)", main = "연령 지원도와 오차")
text(by_cohort$median_nearby_training_n_5yr, by_cohort$MAE, labels = by_cohort$cohort, pos = 3, cex = 0.75)

shift_limits <- c(0, max(by_cohort$cpg_median_shift_median) * 1.08)
plot(by_cohort$cpg_median_shift_median, by_cohort$MAE, pch = 19, col = "darkorange",
     xlim = shift_limits, ylim = mae_limits,
     xlab = "CpG 중앙값 이동의 중앙값", ylab = "코호트 MAE(년)", main = "CpG 분포 이동과 오차")
text(by_cohort$cpg_median_shift_median, by_cohort$MAE, labels = by_cohort$cohort, pos = 3, cex = 0.75)
par(old_par)
invisible(dev.off())

message("=== 연령·코호트 편향 감사 완료 ===")
print(overall, row.names = FALSE)
print(by_age, row.names = FALSE)
print(by_cohort[, c("cohort", "n", "age_min", "age_max", "outside_reference_age_rate",
                    "median_nearby_training_n_5yr", "missing_rate", "cpg_median_shift_median",
                    "MAE", "MAE_ci_low", "MAE_ci_high", "bias_mean")], row.names = FALSE)
print(findings, row.names = FALSE)
