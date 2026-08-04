#!/usr/bin/env Rscript

# 중첩 교차검증에서 선정한 age_density 후보를 8개 학습 코호트만으로 고정한다.
# 외부 검증 데이터 경로를 정의하지 않으며 --confirm-lock 없이는 모델을 저장하지 않는다.

suppressPackageStartupMessages(library(glmnet))
set.seed(20260804)

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
confirm_lock <- "--confirm-lock" %in% args
force <- "--force" %in% args
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

train_path <- file.path(project_dir, "data", "processed", "gse207605_training_bundle.rds")
recommendation_path <- file.path(project_dir, "outputs", "age_bias_mitigation_recommendation.csv")
decision_path <- file.path(project_dir, "outputs", "age_bias_mitigation_candidate_decision.csv")
overall_path <- file.path(project_dir, "outputs", "age_bias_mitigation_metrics_overall.csv")
output_dir <- file.path(project_dir, "outputs")
model_path <- file.path(project_dir, "data", "processed", "locked_age_density_candidate.rds")
tuning_path <- file.path(output_dir, "age_density_candidate_tuning.csv")
selection_path <- file.path(output_dir, "age_density_candidate_selection.csv")
coefficient_path <- file.path(output_dir, "age_density_candidate_nonzero_coefficients.csv")
manifest_path <- file.path(output_dir, "age_density_candidate_lock_manifest.csv")
checkpoint_path <- file.path(output_dir, "age_density_candidate_tuning_checkpoint.rds")

required_paths <- c(train_path, recommendation_path, decision_path, overall_path)
for (path in required_paths) if (!file.exists(path)) stop("필수 입력 파일이 없습니다: ", path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

train <- readRDS(train_path)
x <- train$x
y <- as.numeric(train$y)
group <- as.character(train$group)
sample_id <- as.character(train$sample_id)
features <- as.character(train$feature_names)

if (!is.matrix(x) || nrow(x) != length(y) || length(y) != length(group) || length(y) != length(sample_id)) {
  stop("학습 bundle의 행 수와 메타데이터 길이가 일치하지 않습니다.")
}
if (nrow(x) != 5541L || ncol(x) != 869L) stop("고정 학습 bundle은 5,541명 x 869 CpG여야 합니다.")
if (!identical(colnames(x), features)) stop("CpG 이름 또는 순서가 feature_names와 다릅니다.")
if (anyDuplicated(sample_id) || anyDuplicated(features)) stop("sample_id 또는 CpG 이름에 중복이 있습니다.")
if (any(!is.finite(y)) || any(x < 0 | x > 1, na.rm = TRUE)) stop("연령 또는 beta 값 범위가 올바르지 않습니다.")
cohorts <- sort(unique(group))
if (length(cohorts) != 8L) stop("고정 학습 코호트는 8개여야 합니다.")

recommendation <- read.csv(recommendation_path, stringsAsFactors = FALSE, check.names = FALSE)
decision <- read.csv(decision_path, stringsAsFactors = FALSE, check.names = FALSE)
overall <- read.csv(overall_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(recommendation) != 1L || recommendation$selected_variant[[1]] != "age_density") {
  stop("사전 비교 실험의 선택 후보가 age_density가 아닙니다.")
}
if (!"external_validation_used" %in% names(recommendation) ||
    toupper(as.character(recommendation$external_validation_used[[1]])) != "FALSE") {
  stop("후보 추천 파일의 외부 검증 미사용 표기가 올바르지 않습니다.")
}
if (!"baseline_reproduction_max_abs_diff" %in% names(recommendation) ||
    !is.finite(recommendation$baseline_reproduction_max_abs_diff[[1]]) ||
    recommendation$baseline_reproduction_max_abs_diff[[1]] > 1e-8) {
  stop("기준 모델 OOF 예측 재현 검증을 통과하지 못했습니다.")
}
selected_decision <- decision[decision$candidate == "age_density", , drop = FALSE]
if (nrow(selected_decision) != 1L || !isTRUE(as.logical(selected_decision$adopt[[1]]))) {
  stop("age_density 후보가 사전 채택 조건을 통과하지 않았습니다.")
}
selected_overall <- overall[overall$variant == "age_density", , drop = FALSE]
if (nrow(selected_overall) != 1L || !is.finite(selected_overall$MAE[[1]])) stop("age_density OOF 지표가 없습니다.")

age_breaks <- c(-Inf, 20, 40, 60, 80, Inf)
age_labels <- c("<20", "20-39", "40-59", "60-79", ">=80")
weight_cap <- 5
alphas <- c(0.05, 0.1, 0.5)
# outcome과 무관하게 미리 고정한 lambda 격자다.
lambdas <- 10^seq(2, -4, length.out = 80L)

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

normalize_capped_weights <- function(raw_weights, cap) {
  target_sum <- length(raw_weights)
  weights <- rep(NA_real_, length(raw_weights))
  capped <- rep(FALSE, length(raw_weights))
  repeat {
    uncapped <- !capped
    remaining_sum <- target_sum - sum(weights[capped], na.rm = TRUE)
    scale <- remaining_sum / sum(raw_weights[uncapped])
    proposed <- raw_weights[uncapped] * scale
    exceeds <- proposed > cap
    if (!any(exceeds)) {
      weights[uncapped] <- proposed
      break
    }
    new_capped <- which(uncapped)[exceeds]
    weights[new_capped] <- cap
    capped[new_capped] <- TRUE
  }
  if (abs(sum(weights) - target_sum) > 1e-8 || max(weights) > cap + 1e-10) {
    stop("연령 밀도 가중치 정규화에 실패했습니다.")
  }
  weights
}

make_age_density_weights <- function(ages) {
  bins <- cut(ages, breaks = age_breaks, labels = age_labels, right = FALSE)
  counts <- table(bins)
  reference <- median(as.numeric(counts[counts > 0]))
  raw <- reference / as.numeric(counts[as.character(bins)])
  normalize_capped_weights(raw, weight_cap)
}

full_weights <- make_age_density_weights(y)
if (dry_run) {
  cat("후보: age_density\n")
  cat("학습 표본:", nrow(x), "/ CpG:", ncol(x), "/ 코호트:", length(cohorts), "\n")
  cat("예정 적합 횟수:", length(alphas) * length(cohorts) + 1L, "\n")
  cat("가중치 최소/중앙/최대/합:", min(full_weights), median(full_weights), max(full_weights), sum(full_weights), "\n")
  cat("외부 검증 데이터 경로: 정의되지 않음\n")
  cat("dry-run: 모델과 산출물을 저장하지 않았습니다.\n")
  quit(status = 0)
}

if (!confirm_lock) {
  stop("후보 모델을 저장하려면 --confirm-lock을 명시하세요. 외부 데이터는 사용되지 않습니다.")
}
if (!force && file.exists(model_path)) {
  stop("잠금 후보 모델이 이미 있습니다: ", model_path,
       "\n재생성 근거가 있을 때만 --force를 사용하세요.")
}

signature <- list(
  training_bundle_md5 = unname(tools::md5sum(train_path)),
  sample_count = nrow(x), feature_count = ncol(x), cohorts = cohorts,
  alphas = alphas, lambdas = lambdas, weighting_scheme = "age_density", weight_cap = weight_cap
)
checkpoint <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else list(signature = signature, by_alpha = list())
if (!identical(checkpoint$signature, signature)) stop("기존 tuning checkpoint의 설정이 현재 실행과 다릅니다.")

for (alpha in alphas) {
  key <- format(alpha, scientific = FALSE)
  if (key %in% names(checkpoint$by_alpha)) {
    message("checkpoint 재사용: alpha=", key)
    next
  }
  prediction_grid <- matrix(NA_real_, nrow = length(y), ncol = length(lambdas))
  for (held_out in cohorts) {
    validation <- which(group == held_out)
    fitting <- which(group != held_out)
    medians <- column_medians(x[fitting, , drop = FALSE])
    x_fit <- impute_with(x[fitting, , drop = FALSE], medians)
    x_validation <- impute_with(x[validation, , drop = FALSE], medians)
    fit <- glmnet(
      x_fit, y[fitting], alpha = alpha, lambda = lambdas,
      weights = make_age_density_weights(y[fitting]), standardize = TRUE
    )
    prediction_grid[validation, ] <- predict(fit, newx = x_validation, s = lambdas)
  }
  if (anyNA(prediction_grid)) stop("LOCO 예측에 결측값이 있습니다: alpha=", alpha)
  mae <- colMeans(abs(prediction_grid - y))
  checkpoint$by_alpha[[key]] <- data.frame(
    alpha = alpha, lambda = lambdas, loco_tuning_mae = mae, stringsAsFactors = FALSE
  )
  saveRDS(checkpoint, checkpoint_path)
  message(sprintf("alpha=%.2f 완료 / 최저 LOCO tuning MAE=%.4f", alpha, min(mae)))
}

tuning <- do.call(rbind, checkpoint$by_alpha)
row.names(tuning) <- NULL
best <- tuning[which.min(tuning$loco_tuning_mae), , drop = FALSE]

medians <- column_medians(x)
x_fit <- impute_with(x, medians)
final_fit <- glmnet(
  x_fit, y, alpha = best$alpha[[1]], lambda = lambdas,
  weights = full_weights, standardize = TRUE
)
coefficient_matrix <- as.matrix(coef(final_fit, s = best$lambda[[1]]))
coefficient_table <- data.frame(
  feature = rownames(coefficient_matrix), coefficient = as.numeric(coefficient_matrix[, 1]),
  stringsAsFactors = FALSE
)
nonzero <- coefficient_table[coefficient_table$coefficient != 0, , drop = FALSE]
nonzero_count <- sum(nonzero$feature != "(Intercept)")

locked_model <- list(
  model = final_fit,
  alpha = best$alpha[[1]], lambda = best$lambda[[1]],
  weighting_scheme = "age_density", weight_cap = weight_cap,
  age_breaks = age_breaks, age_labels = age_labels,
  imputation_medians = medians, feature_names = features,
  training_samples = length(y), training_cohorts = cohorts,
  training_age_min = min(y), training_age_max = max(y),
  training_bundle_md5 = signature$training_bundle_md5,
  recommendation_md5 = unname(tools::md5sum(recommendation_path)),
  decision_md5 = unname(tools::md5sum(decision_path)),
  overall_metrics_md5 = unname(tools::md5sum(overall_path)),
  selection_metric = "training-only leave-one-cohort-out tuning MAE",
  selection_mae = best$loco_tuning_mae[[1]],
  nested_cv_candidate_mae = selected_overall$MAE[[1]],
  nonzero_coefficients = nonzero_count,
  external_validation_used = FALSE,
  lock_version = "age-density-v1"
)

saveRDS(locked_model, model_path)
model_md5 <- unname(tools::md5sum(model_path))
selection <- data.frame(
  candidate = "age_density", training_samples = length(y), training_cohorts = length(cohorts),
  feature_count = length(features), selected_alpha = locked_model$alpha,
  selected_lambda = locked_model$lambda, loco_tuning_mae = locked_model$selection_mae,
  nested_cv_candidate_mae = locked_model$nested_cv_candidate_mae,
  nonzero_coefficients = nonzero_count, weight_cap = weight_cap,
  external_validation_used = FALSE, stringsAsFactors = FALSE
)
manifest <- data.frame(
  lock_version = locked_model$lock_version,
  model_file = "data/processed/locked_age_density_candidate.rds",
  model_md5 = model_md5,
  training_bundle_md5 = signature$training_bundle_md5,
  recommendation_md5 = locked_model$recommendation_md5,
  decision_md5 = locked_model$decision_md5,
  overall_metrics_md5 = locked_model$overall_metrics_md5,
  training_samples = length(y), training_cohorts = paste(cohorts, collapse = ";"),
  training_age_min = min(y), training_age_max = max(y), feature_count = length(features),
  weighting_scheme = "age_density", weight_cap = weight_cap,
  selected_alpha = locked_model$alpha, selected_lambda = locked_model$lambda,
  glmnet_version = as.character(packageVersion("glmnet")),
  external_validation_used = FALSE,
  locked_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)

write.csv(tuning, tuning_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(selection, selection_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(nonzero, coefficient_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(manifest, manifest_path, row.names = FALSE, fileEncoding = "UTF-8")

message("=== age_density 학습 전용 후보 잠금 완료 ===")
print(selection, row.names = FALSE)
print(manifest, row.names = FALSE)
