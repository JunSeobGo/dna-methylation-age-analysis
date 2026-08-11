#!/usr/bin/env Rscript

# 외부 검증셋을 사용하지 않고 코호트 단위 중첩 교차검증으로 연령 편향 개선 후보를 비교한다.
# 비교 후보: 기본 가중치, 연령 밀도 가중치, 코호트·연령 혼합 가중치,
#            inner OOF 예측으로만 학습한 기본 모델 calibration.

suppressPackageStartupMessages(library(glmnet))
set.seed(20260803)

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

bundle_path <- file.path(project_dir, "data", "processed", "gse207605_training_bundle.rds")
baseline_prediction_path <- file.path(project_dir, "outputs", "gse207605_cv_oof_predictions.csv")
output_dir <- file.path(project_dir, "outputs")
overall_path <- file.path(output_dir, "age_bias_mitigation_metrics_overall.csv")
age_path <- file.path(output_dir, "age_bias_mitigation_metrics_by_agebin.csv")
cohort_path <- file.path(output_dir, "age_bias_mitigation_metrics_by_cohort.csv")
fold_path <- file.path(output_dir, "age_bias_mitigation_fold_params.csv")
prediction_path <- file.path(output_dir, "age_bias_mitigation_oof_predictions.csv")
decision_path <- file.path(output_dir, "age_bias_mitigation_candidate_decision.csv")
recommendation_path <- file.path(output_dir, "age_bias_mitigation_recommendation.csv")
bootstrap_path <- file.path(output_dir, "age_bias_mitigation_cluster_bootstrap.csv")
plot_path <- file.path(output_dir, "age_bias_mitigation_diagnostics.png")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(bundle_path)) stop("학습 bundle이 없습니다: ", bundle_path)

bundle <- readRDS(bundle_path)
x <- bundle$x
y <- as.numeric(bundle$y)
group <- as.character(bundle$group)
sample_id <- as.character(bundle$sample_id)
features <- as.character(bundle$feature_names)

if (!is.matrix(x) || nrow(x) != length(y) || length(y) != length(group) || length(y) != length(sample_id)) {
  stop("학습 bundle의 행 수와 메타데이터 길이가 일치하지 않습니다.")
}
if (!identical(colnames(x), features) || ncol(x) != 869L) stop("고정된 869개 CpG 순서가 일치하지 않습니다.")
if (anyDuplicated(sample_id) || anyDuplicated(features)) stop("sample_id 또는 CpG 이름에 중복이 있습니다.")
if (any(!is.finite(y))) stop("학습 연령에 결측값 또는 무한값이 있습니다.")
cohorts <- sort(unique(group))
if (length(cohorts) != 8L) stop("이 실험은 고정된 8개 학습 코호트를 전제로 합니다.")

# 기존 standard·balanced 중첩 CV의 16개 outer fold에서 실제 선택된 alpha만 사용한다.
# 모두 0.05, 0.1, 0.5 중 하나였으며, 외부 검증 정보는 이 축소에 사용하지 않았다.
alphas <- c(0.05, 0.1, 0.5)
weight_schemes <- c("standard", "age_density", "mixed")
variants <- c("standard", "age_density", "mixed", "calibrated_standard")
age_breaks <- c(-Inf, 20, 40, 60, 80, Inf)
age_labels <- c("<20", "20-39", "40-59", "60-79", ">=80")
weight_cap <- 5

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

# 평균 가중치는 1로 유지하면서 최종 개별 가중치가 cap을 넘지 않도록 정규화한다.
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
    stop("가중치 상한 정규화에 실패했습니다.")
  }
  weights
}

make_weights <- function(groups, ages, scheme) {
  if (scheme == "standard") return(rep(1, length(ages)))

  bins <- cut(ages, breaks = age_breaks, labels = age_labels, right = FALSE)
  age_counts <- table(bins)
  nonzero_age_counts <- as.numeric(age_counts[age_counts > 0])
  age_reference <- median(nonzero_age_counts)
  age_weights <- age_reference / as.numeric(age_counts[as.character(bins)])

  if (scheme == "age_density") {
    weights <- age_weights
  } else if (scheme == "mixed") {
    cohort_counts <- table(groups)
    cohort_reference <- median(as.numeric(cohort_counts))
    cohort_weights <- cohort_reference / as.numeric(cohort_counts[groups])
    weights <- sqrt(age_weights * cohort_weights)
  } else {
    stop("지원하지 않는 가중치 방식입니다: ", scheme)
  }
  normalize_capped_weights(weights, weight_cap)
}

metric_values <- function(actual, predicted) {
  residual <- predicted - actual
  denominator <- sum((actual - mean(actual))^2)
  c(
    MAE = mean(abs(residual)), RMSE = sqrt(mean(residual^2)),
    MedAE = median(abs(residual)),
    R2 = if (denominator > 0) 1 - sum(residual^2) / denominator else NA_real_,
    bias_mean = mean(residual),
    residual_age_slope = if (length(unique(actual)) > 1L) as.numeric(coef(lm(residual ~ actual))[2]) else NA_real_
  )
}

if (dry_run) {
  design <- data.frame(
    scheme = weight_schemes,
    approximate_fits = length(cohorts) * (length(alphas) * length(cohorts) + 1L),
    stringsAsFactors = FALSE
  )
  weight_profile <- do.call(rbind, lapply(weight_schemes, function(scheme) {
    weights <- make_weights(group, y, scheme)
    data.frame(scheme = scheme, min_weight = min(weights), median_weight = median(weights),
               max_weight = max(weights), weight_sum = sum(weights))
  }))
  message("dry-run: 입력 구조와 가중치·중첩 교차검증 설계만 점검했습니다. 모델은 학습하지 않습니다.")
  cat("표본:", length(y), "명 / CpG:", ncol(x), "개 / 코호트:", length(cohorts), "개\n")
  cat("외부 검증셋은 입력 경로에 포함되지 않습니다.\n")
  print(design, row.names = FALSE)
  print(weight_profile, row.names = FALSE)
  quit(status = 0)
}

run_scheme <- function(scheme) {
  message("== 가중치 방식: ", scheme, " ==")
  checkpoint_path <- file.path(output_dir, paste0("age_bias_mitigation_checkpoint_", scheme, ".rds"))
  progress_path <- file.path(output_dir, paste0("age_bias_mitigation_progress_", scheme, ".csv"))
  signature <- list(
    sample_count = length(y), feature_count = ncol(x), cohorts = cohorts,
    alphas = alphas, weight_cap = weight_cap, scheme = scheme
  )
  checkpoint <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else list(signature = signature, folds = list())
  if (!identical(checkpoint$signature, signature)) {
    stop("기존 checkpoint의 실험 설정이 현재 코드와 다릅니다: ", checkpoint_path)
  }

  for (outer_number in seq_along(cohorts)) {
    held_out <- cohorts[[outer_number]]
    if (held_out %in% names(checkpoint$folds)) {
      message("  [", scheme, "] checkpoint 재사용: ", held_out)
      next
    }
    test_index <- which(group == held_out)
    train_index <- which(group != held_out)
    x_train <- x[train_index, , drop = FALSE]
    y_train <- y[train_index]
    group_train <- group[train_index]
    x_test <- x[test_index, , drop = FALSE]
    inner_cohorts <- sort(unique(group_train))

    outer_medians <- column_medians(x_train)
    x_train_imputed <- impute_with(x_train, outer_medians)
    x_test_imputed <- impute_with(x_test, outer_medians)
    outer_weights <- make_weights(group_train, y_train, scheme)

    best <- list(mae = Inf, alpha = NA_real_, lambda = NA_real_, inner_oof = NULL)

    for (alpha in alphas) {
      lambda_grid <- glmnet(
        x_train_imputed, y_train, alpha = alpha,
        weights = outer_weights, standardize = TRUE
      )$lambda
      inner_predictions <- matrix(NA_real_, nrow = length(y_train), ncol = length(lambda_grid))

      for (inner_held_out in inner_cohorts) {
        validation <- which(group_train == inner_held_out)
        fitting <- which(group_train != inner_held_out)
        inner_medians <- column_medians(x_train[fitting, , drop = FALSE])
        x_fitting <- impute_with(x_train[fitting, , drop = FALSE], inner_medians)
        x_validation <- impute_with(x_train[validation, , drop = FALSE], inner_medians)
        inner_weights <- make_weights(group_train[fitting], y_train[fitting], scheme)
        fit <- glmnet(
          x_fitting, y_train[fitting], alpha = alpha, lambda = lambda_grid,
          weights = inner_weights, standardize = TRUE
        )
        prediction <- predict(fit, newx = x_validation)
        inner_predictions[validation, seq_len(ncol(prediction))] <- prediction
      }

      valid <- colSums(is.na(inner_predictions)) == 0
      if (!any(valid)) next
      mae_by_lambda <- rep(Inf, length(lambda_grid))
      mae_by_lambda[valid] <- colMeans(abs(inner_predictions[, valid, drop = FALSE] - y_train))
      best_index <- which.min(mae_by_lambda)
      if (mae_by_lambda[[best_index]] < best$mae) {
        best <- list(
          mae = mae_by_lambda[[best_index]], alpha = alpha,
          lambda = lambda_grid[[best_index]],
          inner_oof = inner_predictions[, best_index]
        )
      }
    }

    if (!is.finite(best$mae) || anyNA(best$inner_oof)) stop("유효한 inner OOF 후보를 찾지 못했습니다: ", held_out)
    final_fit <- glmnet(
      x_train_imputed, y_train, alpha = best$alpha,
      weights = outer_weights, standardize = TRUE
    )
    raw_prediction <- as.numeric(predict(final_fit, newx = x_test_imputed, s = best$lambda))
    calibration_intercept <- NA_real_
    calibration_slope <- NA_real_
    calibrated_prediction <- NULL
    if (scheme == "standard") {
      calibration_data <- data.frame(actual = y_train, predicted = best$inner_oof)
      calibration_fit <- lm(actual ~ predicted, data = calibration_data)
      calibration_intercept <- as.numeric(coef(calibration_fit)[1])
      calibration_slope <- as.numeric(coef(calibration_fit)[2])
      calibrated_prediction <- calibration_intercept + calibration_slope * raw_prediction
    }

    nonzero <- sum(as.numeric(coef(final_fit, s = best$lambda))[-1] != 0)
    fold_row <- data.frame(
      scheme = scheme, outer_cohort = held_out, n_test = length(test_index),
      alpha = best$alpha, lambda = best$lambda, inner_cv_mae = best$mae,
      nonzero = nonzero, weight_min = min(outer_weights),
      weight_median = median(outer_weights), weight_max = max(outer_weights),
      calibration_intercept = calibration_intercept,
      calibration_slope = calibration_slope,
      stringsAsFactors = FALSE
    )
    checkpoint$folds[[held_out]] <- list(
      test_index = test_index, raw_prediction = raw_prediction,
      calibrated_prediction = calibrated_prediction, params = fold_row
    )
    saveRDS(checkpoint, checkpoint_path)
    write.csv(
      data.frame(
        scheme = scheme, completed_folds = length(checkpoint$folds), total_folds = length(cohorts),
        last_completed_cohort = held_out,
        status = if (length(checkpoint$folds) == length(cohorts)) "completed" else "running",
        updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), stringsAsFactors = FALSE
      ),
      progress_path, row.names = FALSE, fileEncoding = "UTF-8"
    )
    message(sprintf("  [%s] test=%s n=%d alpha=%.2f lambda=%.4f nonzero=%d",
                    scheme, held_out, length(test_index), best$alpha, best$lambda, nonzero))
  }

  oof_raw <- rep(NA_real_, length(y))
  oof_calibrated <- if (scheme == "standard") rep(NA_real_, length(y)) else NULL
  fold_rows <- list()
  for (held_out in cohorts) {
    fold_result <- checkpoint$folds[[held_out]]
    if (is.null(fold_result)) stop("checkpoint에 outer fold가 없습니다: ", held_out)
    oof_raw[fold_result$test_index] <- fold_result$raw_prediction
    if (scheme == "standard") oof_calibrated[fold_result$test_index] <- fold_result$calibrated_prediction
    fold_rows[[held_out]] <- fold_result$params
  }
  if (anyNA(oof_raw)) stop("outer OOF 예측에 결측값이 있습니다: ", scheme)
  if (scheme == "standard" && anyNA(oof_calibrated)) stop("calibrated OOF 예측에 결측값이 있습니다.")
  list(raw = oof_raw, calibrated = oof_calibrated, params = do.call(rbind, fold_rows))
}

results <- lapply(weight_schemes, run_scheme)
names(results) <- weight_schemes
predictions <- list(
  standard = results$standard$raw,
  age_density = results$age_density$raw,
  mixed = results$mixed$raw,
  calibrated_standard = results$standard$calibrated
)

# 기존 기준 스크립트의 standard 예측과 새 구현이 일치해야 비교 기준이 흔들리지 않는다.
baseline_reproduction_max_abs_diff <- NA_real_
if (file.exists(baseline_prediction_path)) {
  prior <- read.csv(baseline_prediction_path, stringsAsFactors = FALSE)
  if (!all(c("sample_id", "pred_standard") %in% names(prior))) stop("기존 기준 예측 파일 형식이 올바르지 않습니다.")
  prior <- prior[match(sample_id, prior$sample_id), , drop = FALSE]
  if (anyNA(prior$pred_standard)) stop("기존 기준 예측과 sample_id가 일치하지 않습니다.")
  baseline_reproduction_max_abs_diff <- max(abs(prior$pred_standard - predictions$standard))
  if (baseline_reproduction_max_abs_diff > 1e-8) {
    stop("기준 모델 재현 실패: 최대 예측 차이 = ", baseline_reproduction_max_abs_diff)
  }
}

metrics_overall <- function(variant, predicted) {
  data.frame(variant = variant, n = length(y), t(metric_values(y, predicted)), stringsAsFactors = FALSE)
}

metrics_by_age <- function(variant, predicted) {
  bins <- cut(y, breaks = age_breaks, labels = age_labels, right = FALSE)
  do.call(rbind, lapply(age_labels, function(bin) {
    index <- which(bins == bin)
    data.frame(variant = variant, age_bin = bin, n = length(index),
               t(metric_values(y[index], predicted[index])), stringsAsFactors = FALSE)
  }))
}

metrics_by_cohort <- function(variant, predicted) {
  do.call(rbind, lapply(cohorts, function(cohort) {
    index <- which(group == cohort)
    data.frame(variant = variant, cohort = cohort, n = length(index),
               age_min = min(y[index]), age_max = max(y[index]),
               t(metric_values(y[index], predicted[index])), stringsAsFactors = FALSE)
  }))
}

overall <- do.call(rbind, lapply(variants, function(variant) metrics_overall(variant, predictions[[variant]])))
by_age <- do.call(rbind, lapply(variants, function(variant) metrics_by_age(variant, predictions[[variant]])))
by_cohort <- do.call(rbind, lapply(variants, function(variant) metrics_by_cohort(variant, predictions[[variant]])))
fold_params <- do.call(rbind, lapply(weight_schemes, function(scheme) results[[scheme]]$params))

oof <- data.frame(
  sample_id = sample_id, cohort = group, age_true = y,
  pred_standard = predictions$standard,
  pred_age_density = predictions$age_density,
  pred_mixed = predictions$mixed,
  pred_calibrated_standard = predictions$calibrated_standard,
  stringsAsFactors = FALSE
)

baseline_overall <- overall[overall$variant == "standard", ]
baseline_young <- by_age$MAE[by_age$variant == "standard" & by_age$age_bin == "<20"]
baseline_old <- by_age$MAE[by_age$variant == "standard" & by_age$age_bin == ">=80"]
baseline_worst <- max(by_cohort$MAE[by_cohort$variant == "standard"])
candidates <- setdiff(variants, "standard")

decision <- do.call(rbind, lapply(candidates, function(candidate) {
  candidate_overall <- overall[overall$variant == candidate, ]
  candidate_young <- by_age$MAE[by_age$variant == candidate & by_age$age_bin == "<20"]
  candidate_old <- by_age$MAE[by_age$variant == candidate & by_age$age_bin == ">=80"]
  candidate_worst <- max(by_cohort$MAE[by_cohort$variant == candidate])
  overall_delta <- candidate_overall$MAE - baseline_overall$MAE
  young_improvement <- baseline_young - candidate_young
  old_improvement <- baseline_old - candidate_old
  worst_delta <- candidate_worst - baseline_worst
  pass_overall <- overall_delta <= 0.20
  pass_extreme <- max(young_improvement, old_improvement) >= 0.50
  pass_worst <- worst_delta <= 0
  data.frame(
    candidate = candidate,
    overall_MAE = candidate_overall$MAE, overall_MAE_delta = overall_delta,
    under20_MAE = candidate_young, under20_MAE_improvement = young_improvement,
    age80plus_MAE = candidate_old, age80plus_MAE_improvement = old_improvement,
    worst_cohort_MAE = candidate_worst, worst_cohort_MAE_delta = worst_delta,
    pass_overall_guardrail = pass_overall,
    pass_extreme_improvement = pass_extreme,
    pass_worst_cohort_guardrail = pass_worst,
    adopt = pass_overall && pass_extreme && pass_worst,
    stringsAsFactors = FALSE
  )
}))

eligible <- decision[decision$adopt, , drop = FALSE]
if (nrow(eligible)) {
  extreme_mean <- (eligible$under20_MAE + eligible$age80plus_MAE) / 2
  selected_variant <- eligible$candidate[[which.min(extreme_mean)]]
  recommendation_reason <- "사전 채택 기준을 모두 통과한 후보 중 연령 양끝 평균 MAE가 가장 낮음"
} else {
  selected_variant <- "standard"
  recommendation_reason <- "사전 채택 기준을 모두 통과한 후보가 없어 기준 모델 유지"
}
recommendation <- data.frame(
  selected_variant = selected_variant,
  recommendation_reason = recommendation_reason,
  baseline_reproduction_max_abs_diff = baseline_reproduction_max_abs_diff,
  external_validation_used = FALSE,
  stringsAsFactors = FALSE
)

cluster_bootstrap <- function(candidate, iterations = 2000L) {
  set.seed(20260803L + match(candidate, candidates))
  baseline <- predictions$standard
  alternative <- predictions[[candidate]]
  bins <- cut(y, breaks = age_breaks, labels = age_labels, right = FALSE)
  estimates <- matrix(NA_real_, nrow = iterations, ncol = 3L)
  colnames(estimates) <- c("overall", "under20", "age80plus")
  for (iteration in seq_len(iterations)) {
    sampled_cohorts <- sample(cohorts, length(cohorts), replace = TRUE)
    index <- unlist(lapply(sampled_cohorts, function(cohort) which(group == cohort)), use.names = FALSE)
    estimates[iteration, "overall"] <- mean(abs(baseline[index] - y[index])) - mean(abs(alternative[index] - y[index]))
    young_index <- index[bins[index] == "<20"]
    old_index <- index[bins[index] == ">=80"]
    if (length(young_index)) {
      estimates[iteration, "under20"] <- mean(abs(baseline[young_index] - y[young_index])) -
        mean(abs(alternative[young_index] - y[young_index]))
    }
    if (length(old_index)) {
      estimates[iteration, "age80plus"] <- mean(abs(baseline[old_index] - y[old_index])) -
        mean(abs(alternative[old_index] - y[old_index]))
    }
  }
  do.call(rbind, lapply(colnames(estimates), function(segment) {
    values <- estimates[, segment]
    values <- values[is.finite(values)]
    data.frame(
      candidate = candidate, segment = segment, valid_iterations = length(values),
      MAE_improvement_mean = mean(values),
      MAE_improvement_ci_low = unname(quantile(values, 0.025)),
      MAE_improvement_ci_high = unname(quantile(values, 0.975)),
      stringsAsFactors = FALSE
    )
  }))
}
bootstrap <- do.call(rbind, lapply(candidates, cluster_bootstrap))

write.csv(overall, overall_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(by_age, age_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(by_cohort, cohort_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fold_params, fold_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(oof, prediction_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(decision, decision_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(recommendation, recommendation_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bootstrap, bootstrap_path, row.names = FALSE, fileEncoding = "UTF-8")

png(plot_path, width = 1400, height = 1000, res = 130)
old_par <- par(mfrow = c(2, 2), mar = c(5.5, 4.5, 3.2, 1))
variant_labels <- c("기본", "연령", "혼합", "보정")
barplot(overall$MAE, names.arg = variant_labels, col = "steelblue", ylim = c(0, max(overall$MAE) * 1.15),
        ylab = "MAE(년)", main = "전체 중첩 교차검증 MAE")
text(seq_along(variants) * 1.2 - 0.5, overall$MAE, labels = sprintf("%.2f", overall$MAE), pos = 3)

extreme_matrix <- rbind(
  by_age$MAE[match(paste(variants, "<20"), paste(by_age$variant, by_age$age_bin))],
  by_age$MAE[match(paste(variants, ">=80"), paste(by_age$variant, by_age$age_bin))]
)
barplot(extreme_matrix, beside = TRUE, names.arg = variant_labels, col = c("darkorange", "tomato"),
        ylab = "MAE(년)", main = "연령 양끝 MAE", legend.text = c("<20", ">=80"),
        args.legend = list(x = "topright", bty = "n"))

worst_values <- vapply(variants, function(variant) max(by_cohort$MAE[by_cohort$variant == variant]), numeric(1))
barplot(worst_values, names.arg = variant_labels, col = "slateblue", ylim = c(0, max(worst_values) * 1.15),
        ylab = "최악 코호트 MAE(년)", main = "코호트 강건성")

plot(y, predictions$standard - y, pch = 16, cex = 0.3, col = adjustcolor("grey40", alpha.f = 0.2),
     xlab = "실제 나이", ylab = "잔차(예측 - 실제, 년)", main = "후보별 잔차 추세")
abline(h = 0, col = "black", lty = 2)
line_colors <- c("black", "darkorange", "steelblue", "purple")
for (index in seq_along(variants)) {
  lines(lowess(y, predictions[[variants[[index]]]] - y), col = line_colors[[index]], lwd = 2)
}
legend("bottomleft", legend = variant_labels, col = line_colors, lwd = 2, bty = "n")
par(old_par)
invisible(dev.off())

message("=== 연령 편향 개선 후보 비교 완료 ===")
print(overall, row.names = FALSE)
print(by_age, row.names = FALSE)
print(decision, row.names = FALSE)
print(recommendation, row.names = FALSE)
print(bootstrap, row.names = FALSE)
