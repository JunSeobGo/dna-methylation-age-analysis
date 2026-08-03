#!/usr/bin/env Rscript

# 외부 검증셋을 사용하지 않고 8개 학습 코호트 내부에서 표본 수별 학습곡선을 측정한다.
# 각 평가 코호트를 통째로 제외한 뒤, 나머지 코호트에서 연령대별로 층화 표본추출한다.

suppressPackageStartupMessages(library(glmnet))

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

bundle_path <- file.path(project_dir, "data", "processed", "gse207605_training_bundle.rds")
model_path <- file.path(project_dir, "data", "processed", "locked_elastic_net_model.rds")
output_dir <- file.path(project_dir, "outputs")
prediction_path <- file.path(output_dir, "training_learning_curve_predictions.csv")
metric_path <- file.path(output_dir, "training_learning_curve_metrics.csv")
summary_path <- file.path(output_dir, "training_learning_curve_summary.csv")
age_path <- file.path(output_dir, "training_learning_curve_by_agebin.csv")
cohort_path <- file.path(output_dir, "training_learning_curve_by_cohort.csv")
decision_path <- file.path(output_dir, "training_learning_curve_decision.csv")
plot_path <- file.path(output_dir, "training_learning_curve.png")

if (!file.exists(bundle_path)) stop("학습 bundle이 없습니다: ", bundle_path)
if (!file.exists(model_path)) stop("고정 모델이 없습니다: ", model_path)

bundle <- readRDS(bundle_path)
locked <- readRDS(model_path)
x <- bundle$x
y <- as.numeric(bundle$y)
group <- as.character(bundle$group)
sample_id <- as.character(bundle$sample_id)
features <- as.character(bundle$feature_names)

if (!is.matrix(x) || nrow(x) != length(y) || length(y) != length(group)) {
  stop("학습 bundle의 행 수와 결과·코호트 벡터 길이가 일치하지 않습니다.")
}
if (any(!is.finite(y))) stop("학습 연령에 결측값 또는 무한값이 있습니다.")
if (anyDuplicated(sample_id)) stop("학습 bundle에 중복 sample_id가 있습니다.")
if (anyDuplicated(features)) stop("학습 bundle에 중복 CpG가 있습니다.")
if (!identical(colnames(x), features)) stop("학습 bundle의 CpG 이름과 행렬 열 순서가 다릅니다.")
if (!identical(features, as.character(locked$feature_names))) stop("고정 모델과 학습 bundle의 CpG 순서가 다릅니다.")
cohorts <- sort(unique(group))
if (length(cohorts) != 8L) stop("학습곡선은 고정된 8개 학습 코호트를 전제로 합니다.")
if (!identical(sort(as.character(locked$training_cohorts)), cohorts)) {
  stop("고정 모델의 학습 코호트와 현재 bundle의 코호트가 다릅니다.")
}
if (!locked$weighting_scheme %in% c("standard", "balanced")) {
  stop("지원하지 않는 고정 모델 가중치 방식입니다: ", locked$weighting_scheme)
}

fractions <- c(0.25, 0.50, 0.75, 1.00)
repeats_by_fraction <- c(`0.25` = 5L, `0.5` = 5L, `0.75` = 5L, `1` = 1L)
base_seed <- 20260803L
age_levels <- c("<20", "20-39", "40-59", "60-79", ">=80")

age_bin <- function(values) {
  cut(values, breaks = c(-Inf, 20, 40, 60, 80, Inf), labels = age_levels, right = FALSE)
}

column_medians <- function(mat) {
  values <- apply(mat, 2, median, na.rm = TRUE)
  values[!is.finite(values)] <- 0.5
  values
}

impute_with <- function(mat, medians) {
  positions <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(positions)) mat[positions] <- medians[positions[, "col"]]
  mat
}

make_weights <- function(groups, scheme) {
  if (scheme == "standard") return(rep(1, length(groups)))
  sizes <- table(groups)
  weights <- as.numeric(1 / sizes[groups])
  weights * length(groups) / sum(weights)
}

sample_stratified <- function(indices, fraction, seed) {
  if (fraction == 1) return(indices)
  set.seed(seed)
  strata <- interaction(group[indices], age_bin(y[indices]), drop = TRUE)
  cells <- split(indices, strata)
  selected <- unlist(lapply(cells, function(cell) {
    target <- max(1L, floor(length(cell) * fraction))
    sample(cell, size = target, replace = FALSE)
  }), use.names = FALSE)
  sort(selected)
}

metric_values <- function(actual, predicted) {
  residual <- predicted - actual
  denominator <- sum((actual - mean(actual))^2)
  c(
    MAE = mean(abs(residual)),
    RMSE = sqrt(mean(residual^2)),
    MedAE = median(abs(residual)),
    R2 = if (denominator > 0) 1 - sum(residual^2) / denominator else NA_real_,
    bias_mean = mean(residual),
    residual_age_slope = if (length(unique(actual)) > 1L) as.numeric(coef(lm(residual ~ actual))[2]) else NA_real_
  )
}

if (dry_run) {
  design <- expand.grid(fraction = fractions, held_out_cohort = cohorts, stringsAsFactors = FALSE)
  design$repeats <- repeats_by_fraction[as.character(design$fraction)]
  design$planned_fits <- design$repeats
  message("dry-run: 입력 구조와 학습곡선 설계만 점검했습니다. 모델은 학습하지 않습니다.")
  cat("학습 표본:", length(y), "명 / CpG:", ncol(x), "개 / 코호트:", length(cohorts), "개\n")
  cat("고정 설정: alpha =", locked$alpha, ", lambda =", locked$lambda,
      ", weighting =", locked$weighting_scheme, "\n")
  cat("총 적합 횟수:", sum(design$planned_fits), "회\n")
  print(aggregate(planned_fits ~ fraction, design, sum), row.names = FALSE)
  quit(status = 0)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
prediction_rows <- list()
row_number <- 1L

for (fraction in fractions) {
  repeat_count <- repeats_by_fraction[[as.character(fraction)]]
  for (repeat_id in seq_len(repeat_count)) {
    for (cohort_number in seq_along(cohorts)) {
      held_out <- cohorts[[cohort_number]]
      validation <- which(group == held_out)
      available <- which(group != held_out)
      seed <- base_seed + as.integer(fraction * 1000) + repeat_id * 100L + cohort_number
      fitting <- sample_stratified(available, fraction, seed)

      medians <- column_medians(x[fitting, , drop = FALSE])
      x_fit <- impute_with(x[fitting, , drop = FALSE], medians)
      x_validation <- impute_with(x[validation, , drop = FALSE], medians)
      fit <- glmnet(x_fit, y[fitting], alpha = locked$alpha, lambda = locked$lambda,
                    weights = make_weights(group[fitting], locked$weighting_scheme), standardize = TRUE)
      predicted <- as.numeric(predict(fit, newx = x_validation, s = locked$lambda))

      prediction_rows[[row_number]] <- data.frame(
        fraction = fraction, repeat_id = repeat_id, held_out_cohort = held_out,
        training_n = length(fitting), sample_id = sample_id[validation],
        age_true = y[validation], age_predicted = predicted,
        residual = predicted - y[validation], age_bin = as.character(age_bin(y[validation])),
        stringsAsFactors = FALSE
      )
      row_number <- row_number + 1L
      message(sprintf("완료: 비율 %3d%% / 반복 %d / 홀드아웃 %s / 학습 n=%d",
                      round(fraction * 100), repeat_id, held_out, length(fitting)))
    }
  }
}

predictions <- do.call(rbind, prediction_rows)
expected_rows <- length(y) * sum(repeats_by_fraction)
if (nrow(predictions) != expected_rows) stop("예측 행 수가 설계값과 다릅니다.")
prediction_key <- paste(predictions$fraction, predictions$repeat_id, predictions$sample_id, sep = "|")
if (anyDuplicated(prediction_key)) stop("비율·반복·sample_id 조합에 중복 예측이 있습니다.")

metric_groups <- split(predictions, interaction(predictions$fraction, predictions$repeat_id, drop = TRUE))
metrics <- do.call(rbind, lapply(metric_groups, function(part) {
  values <- metric_values(part$age_true, part$age_predicted)
  old <- part[part$age_true >= 80, , drop = FALSE]
  old_values <- if (nrow(old)) metric_values(old$age_true, old$age_predicted) else rep(NA_real_, 6L)
  cohort_mae <- aggregate(abs(residual) ~ held_out_cohort, part, mean)
  data.frame(
    fraction = part$fraction[[1]], repeat_id = part$repeat_id[[1]],
    mean_training_n = mean(tapply(part$training_n, part$held_out_cohort, unique)),
    n_predictions = nrow(part), t(values), n_age_80_plus = nrow(old),
    MAE_age_80_plus = old_values[["MAE"]], bias_age_80_plus = old_values[["bias_mean"]],
    worst_cohort_MAE = max(cohort_mae$`abs(residual)`), stringsAsFactors = FALSE
  )
}))
metrics <- metrics[order(metrics$fraction, metrics$repeat_id), ]

summarise_metric <- function(values) {
  c(mean = mean(values, na.rm = TRUE), sd = if (sum(is.finite(values)) > 1L) sd(values, na.rm = TRUE) else 0)
}
summary_rows <- lapply(fractions, function(fraction) {
  part <- metrics[metrics$fraction == fraction, , drop = FALSE]
  mae <- summarise_metric(part$MAE); rmse <- summarise_metric(part$RMSE)
  old_mae <- summarise_metric(part$MAE_age_80_plus); worst <- summarise_metric(part$worst_cohort_MAE)
  data.frame(
    fraction = fraction, repeats = nrow(part), mean_training_n = mean(part$mean_training_n),
    MAE_mean = mae[["mean"]], MAE_sd = mae[["sd"]], RMSE_mean = rmse[["mean"]], RMSE_sd = rmse[["sd"]],
    MAE_age_80_plus_mean = old_mae[["mean"]], MAE_age_80_plus_sd = old_mae[["sd"]],
    worst_cohort_MAE_mean = worst[["mean"]], worst_cohort_MAE_sd = worst[["sd"]]
  )
})
curve_summary <- do.call(rbind, summary_rows)

age_groups <- split(predictions, interaction(predictions$fraction, predictions$repeat_id,
                                              predictions$age_bin, drop = TRUE))
by_age <- do.call(rbind, lapply(age_groups, function(part) {
  values <- metric_values(part$age_true, part$age_predicted)
  data.frame(fraction = part$fraction[[1]], repeat_id = part$repeat_id[[1]],
             age_bin = part$age_bin[[1]], n = nrow(part), t(values), stringsAsFactors = FALSE)
}))
by_age <- by_age[order(by_age$fraction, by_age$repeat_id, match(by_age$age_bin, age_levels)), ]

cohort_groups <- split(predictions, interaction(predictions$fraction, predictions$repeat_id,
                                                 predictions$held_out_cohort, drop = TRUE))
by_cohort <- do.call(rbind, lapply(cohort_groups, function(part) {
  values <- metric_values(part$age_true, part$age_predicted)
  data.frame(fraction = part$fraction[[1]], repeat_id = part$repeat_id[[1]],
             held_out_cohort = part$held_out_cohort[[1]], n = nrow(part),
             training_n = unique(part$training_n), t(values), stringsAsFactors = FALSE)
}))
by_cohort <- by_cohort[order(by_cohort$fraction, by_cohort$repeat_id, by_cohort$held_out_cohort), ]

mae_75 <- curve_summary$MAE_mean[curve_summary$fraction == 0.75]
mae_100 <- curve_summary$MAE_mean[curve_summary$fraction == 1.00]
old_75 <- curve_summary$MAE_age_80_plus_mean[curve_summary$fraction == 0.75]
old_100 <- curve_summary$MAE_age_80_plus_mean[curve_summary$fraction == 1.00]
plateau_threshold <- 0.10
improvement <- mae_75 - mae_100
decision <- data.frame(
  comparison = "75% 대비 100% 학습", overall_MAE_improvement_years = improvement,
  age_80_plus_MAE_improvement_years = old_75 - old_100,
  practical_plateau_threshold_years = plateau_threshold,
  interpretation = if (improvement < plateau_threshold) {
    "전체 MAE 기준으로 실용적 포화 구간에 가깝습니다. 데이터 추가보다 코호트 편향과 고령층 보강을 우선 검토합니다."
  } else {
    "전체 MAE가 아직 의미 있게 개선됩니다. 동일 조건의 독립 코호트 추가가 성능 향상에 도움 될 가능성이 있습니다."
  }, stringsAsFactors = FALSE
)

write.csv(predictions, prediction_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(metrics, metric_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(curve_summary, summary_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(by_age, age_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(by_cohort, cohort_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(decision, decision_path, row.names = FALSE, fileEncoding = "UTF-8")

png(plot_path, width = 1500, height = 500, res = 130)
old_par <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
panels <- list(
  list(mean = curve_summary$MAE_mean, sd = curve_summary$MAE_sd, ylab = "전체 MAE (년)", title = "전체 표본"),
  list(mean = curve_summary$MAE_age_80_plus_mean, sd = curve_summary$MAE_age_80_plus_sd,
       ylab = "80세 이상 MAE (년)", title = "고령층"),
  list(mean = curve_summary$worst_cohort_MAE_mean, sd = curve_summary$worst_cohort_MAE_sd,
       ylab = "최악 코호트 MAE (년)", title = "코호트 강건성")
)
for (panel in panels) {
  lower <- pmax(0, panel$mean - panel$sd); upper <- panel$mean + panel$sd
  plot(curve_summary$mean_training_n, panel$mean, type = "b", pch = 19, lwd = 2,
       ylim = range(c(lower, upper)), xlab = "평균 학습 표본 수", ylab = panel$ylab, main = panel$title)
  variable <- which(panel$sd > 0)
  if (length(variable)) {
    arrows(curve_summary$mean_training_n[variable], lower[variable],
           curve_summary$mean_training_n[variable], upper[variable],
           angle = 90, code = 3, length = 0.05)
  }
  text(curve_summary$mean_training_n, panel$mean,
       labels = paste0(round(curve_summary$fraction * 100), "%"), pos = 3, cex = 0.8)
}
par(old_par)
invisible(dev.off())

message("=== 학습곡선 측정 완료 ===")
print(curve_summary, row.names = FALSE)
print(decision, row.names = FALSE)
