#!/usr/bin/env Rscript

# 학습 코호트만으로 최종 Elastic Net 설정을 고정한 뒤 잠금 외부 세트 GSE87571을 평가한다.
# --confirm-external-evaluation 없이는 외부 데이터를 읽지 않으며, 결과가 있으면 재실행을 막는다.

suppressPackageStartupMessages(library(glmnet))
set.seed(20260803)

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
confirm_external <- "--confirm-external-evaluation" %in% args
force <- "--force" %in% args
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

train_path <- file.path(project_dir, "data", "processed", "gse207605_training_bundle.rds")
test_path <- file.path(project_dir, "data", "processed", "gse87571_external_validation_bundle.rds")
benchmark_path <- file.path(project_dir, "outputs", "gse207605_cv_metrics_overall.csv")
model_path <- file.path(project_dir, "data", "processed", "locked_elastic_net_model.rds")
selection_path <- file.path(project_dir, "outputs", "final_model_selection.csv")
metric_path <- file.path(project_dir, "outputs", "gse87571_external_metrics_overall.csv")
age_metric_path <- file.path(project_dir, "outputs", "gse87571_external_metrics_by_agebin.csv")
prediction_path <- file.path(project_dir, "outputs", "gse87571_external_predictions.csv")
plot_path <- file.path(project_dir, "outputs", "gse87571_external_diagnostics.png")

for (path in c(train_path, benchmark_path)) if (!file.exists(path)) stop("필수 입력 파일이 없습니다: ", path)
if (!force && file.exists(metric_path)) {
  stop("외부 평가 결과가 이미 있습니다: ", metric_path,
       "\n잠금 외부 검증을 반복하지 않습니다. 재실행 근거가 있을 때만 --force를 사용하세요.")
}

train <- readRDS(train_path)
x <- train$x; y <- train$y; group <- train$group; sample_id <- train$sample_id; features <- train$feature_names
if (!is.matrix(x) || nrow(x) != length(y) || length(y) != length(group)) stop("학습 bundle 구조가 일치하지 않습니다.")
if (ncol(x) != 869L || !identical(colnames(x), features)) stop("학습 feature는 순서가 고정된 869개 CpG여야 합니다.")
if (anyDuplicated(sample_id) || anyDuplicated(features)) stop("학습 bundle에 중복 ID가 있습니다.")
cohorts <- sort(unique(group))
if (length(cohorts) != 8L) stop("8개 학습 코호트가 필요합니다.")

benchmark <- read.csv(benchmark_path, stringsAsFactors = FALSE)
if (!all(c("scheme", "MAE") %in% names(benchmark))) stop("내부 검증 지표 파일 형식이 올바르지 않습니다.")
benchmark <- benchmark[is.finite(benchmark$MAE), , drop = FALSE]
scheme <- benchmark$scheme[[which.min(benchmark$MAE)]]
if (!scheme %in% c("standard", "balanced")) stop("지원하지 않는 가중치 방식입니다: ", scheme)

column_medians <- function(mat) {
  values <- apply(mat, 2, median, na.rm = TRUE)
  values[is.na(values)] <- 0.5
  values
}
impute_with <- function(mat, medians) {
  positions <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(positions)) mat[positions] <- medians[positions[, "col"]]
  mat
}
make_weights <- function(groups, selected_scheme) {
  if (selected_scheme == "standard") return(rep(1, length(groups)))
  size <- table(groups)
  weights <- as.numeric(1 / size[as.character(groups)])
  weights * length(groups) / sum(weights)
}

# outcome과 beta 값으로 자동 생성하지 않는 고정 격자다.
alphas <- c(0.05, 0.1, 0.25, 0.5, 0.75, 1.0)
lambdas <- 10^seq(2, -4, length.out = 80L)

# 모델 선택은 학습 데이터 안에서만 leave-one-cohort-out으로 수행한다.
tuning_rows <- list(); row_number <- 1L
for (alpha in alphas) {
  prediction_grid <- matrix(NA_real_, nrow = length(y), ncol = length(lambdas))
  for (held_out in cohorts) {
    validation <- which(group == held_out)
    fitting <- which(group != held_out)
    fold_medians <- column_medians(x[fitting, , drop = FALSE])
    x_fit <- impute_with(x[fitting, , drop = FALSE], fold_medians)
    x_validation <- impute_with(x[validation, , drop = FALSE], fold_medians)
    fit <- glmnet(x_fit, y[fitting], alpha = alpha, lambda = lambdas,
                  weights = make_weights(group[fitting], scheme), standardize = TRUE)
    prediction_grid[validation, ] <- predict(fit, newx = x_validation, s = lambdas)
  }
  mae <- colMeans(abs(prediction_grid - y))
  for (index in seq_along(lambdas)) {
    tuning_rows[[row_number]] <- data.frame(alpha = alpha, lambda = lambdas[[index]], loco_mae = mae[[index]])
    row_number <- row_number + 1L
  }
}
tuning <- do.call(rbind, tuning_rows)
best <- tuning[which.min(tuning$loco_mae), , drop = FALSE]

if (dry_run) {
  message("dry-run: 외부 bundle을 읽거나 파일을 저장하지 않았습니다.")
  print(best, row.names = FALSE)
  quit(status = 0)
}

# 여기까지는 학습 데이터만 사용한다. 결측 대치값과 하이퍼파라미터를 먼저 저장한다.
medians <- column_medians(x)
x_fit <- impute_with(x, medians)
final_fit <- glmnet(x_fit, y, alpha = best$alpha[[1]], lambda = lambdas,
                    weights = make_weights(group, scheme), standardize = TRUE)
nonzero <- sum(as.numeric(coef(final_fit, s = best$lambda[[1]]))[-1] != 0)
locked_model <- list(
  model = final_fit, alpha = best$alpha[[1]], lambda = best$lambda[[1]],
  weighting_scheme = scheme, imputation_medians = medians, feature_names = features,
  training_samples = length(y), training_cohorts = cohorts,
  selection_metric = "leave-one-cohort-out MAE", selection_mae = best$loco_mae[[1]],
  nonzero_coefficients = nonzero, created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(locked_model, model_path)
selection <- data.frame(
  weighting_scheme = scheme, training_samples = length(y), training_cohorts = length(cohorts),
  selected_alpha = locked_model$alpha, selected_lambda = locked_model$lambda,
  selected_loco_mae = locked_model$selection_mae, nonzero_coefficients = nonzero,
  benchmark_mae = benchmark$MAE[[match(scheme, benchmark$scheme)]]
)
write.csv(selection, selection_path, row.names = FALSE, fileEncoding = "UTF-8")
message("학습 데이터만으로 최종 모델을 고정했습니다: ", model_path)
print(selection, row.names = FALSE)

if (!confirm_external) {
  message("외부 데이터를 읽지 않았습니다. 실제 평가는 --confirm-external-evaluation을 사용하세요.")
  quit(status = 0)
}

# 잠금 외부 평가: 모델 고정 파일을 만든 다음에만 실행한다.
if (!file.exists(test_path)) stop("외부 검증 bundle이 없습니다: ", test_path)
test <- readRDS(test_path)
if (!is.matrix(test$x) || !identical(colnames(test$x), locked_model$feature_names)) stop("외부 CpG 순서가 학습 모델과 다릅니다.")
if (length(intersect(test$sample_id, sample_id))) stop("학습·외부 GSM이 중복됩니다.")
valid <- is.finite(test$y)
if (!any(valid)) stop("외부 연령 라벨이 없습니다.")
x_test <- impute_with(test$x[valid, , drop = FALSE], locked_model$imputation_medians)
y_test <- test$y[valid]
ids_test <- test$sample_id[valid]
predicted <- as.numeric(predict(locked_model$model, newx = x_test, s = locked_model$lambda))
residual <- predicted - y_test
calibration <- lm(y_test ~ predicted)

overall <- data.frame(
  dataset = "GSE87571", n = length(y_test), MAE = mean(abs(residual)),
  RMSE = sqrt(mean(residual^2)), MedAE = median(abs(residual)),
  R2 = 1 - sum(residual^2) / sum((y_test - mean(y_test))^2), bias_mean = mean(residual),
  residual_age_slope = as.numeric(coef(lm(residual ~ y_test))[2]),
  calibration_intercept = as.numeric(coef(calibration)[1]), calibration_slope = as.numeric(coef(calibration)[2]),
  selected_alpha = locked_model$alpha, selected_lambda = locked_model$lambda,
  weighting_scheme = locked_model$weighting_scheme, nonzero_coefficients = nonzero
)
bins <- cut(y_test, breaks = c(-Inf, 20, 40, 60, 80, Inf),
            labels = c("<20", "20-39", "40-59", "60-79", ">=80"), right = FALSE)
by_age <- do.call(rbind, lapply(levels(bins), function(bin) {
  index <- which(bins == bin); error <- residual[index]
  data.frame(dataset = "GSE87571", age_bin = bin, n = length(index), MAE = mean(abs(error)),
             RMSE = sqrt(mean(error^2)), MedAE = median(abs(error)), bias_mean = mean(error))
}))
predictions <- data.frame(sample_id = ids_test, cohort = "GSE87571", age_true = y_test,
                          age_predicted = predicted, residual = residual)
write.csv(overall, metric_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(by_age, age_metric_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(predictions, prediction_path, row.names = FALSE, fileEncoding = "UTF-8")

png(plot_path, width = 1200, height = 600, res = 120)
old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(y_test, predicted, pch = 16, cex = 0.55, col = rgb(0, 0, 0, 0.3),
     xlab = "실제 나이", ylab = "예측 나이", main = sprintf("GSE87571 외부 검증 (MAE %.2f)", overall$MAE))
abline(0, 1, col = "red", lwd = 2)
plot(y_test, residual, pch = 16, cex = 0.55, col = rgb(0, 0, 0, 0.3),
     xlab = "실제 나이", ylab = "잔차 (예측 - 실제)", main = "GSE87571 외부 검증 잔차")
abline(h = 0, col = "red", lwd = 2); lines(lowess(y_test, residual), col = "blue", lwd = 2)
par(old_par); invisible(dev.off())

message("=== GSE87571 잠금 외부 검증 결과 ===")
print(overall, row.names = FALSE)
print(by_age, row.names = FALSE)
