#!/usr/bin/env Rscript

# SATSA 품질 게이트를 통과한 동일 표본에 두 잠금 모델을 한 번 적용한다.
# dry-run은 외부 bundle을 읽지 않으며, 실제 평가는 명시적 확인 인자가 있어야 실행된다.

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
confirm_external <- "--confirm-external-evaluation" %in% args
force <- "--force" %in% args
reason_arg <- grep("^--rerun-reason=", args, value = TRUE)
rerun_reason <- if (length(reason_arg)) sub("^--rerun-reason=", "", reason_arg[[1]]) else ""
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

expected_samples <- 447L
expected_features <- 869L
expected_age80plus <- 81L

policy_path <- file.path(project_dir, "config", "satsa_evaluation_policy.csv")
bundle_path <- file.path(project_dir, "data", "processed", "satsa_external_validation_bundle.rds")
candidate_model_path <- file.path(project_dir, "data", "processed", "locked_age_density_candidate.rds")
standard_model_path <- file.path(project_dir, "data", "processed", "locked_elastic_net_model.rds")
candidate_lock_path <- file.path(project_dir, "outputs", "age_density_candidate_lock_manifest.csv")
output_dir <- file.path(project_dir, "outputs")
metric_path <- file.path(output_dir, "satsa_locked_metrics.csv")
prediction_path <- file.path(output_dir, "satsa_locked_predictions.csv")
bootstrap_path <- file.path(output_dir, "satsa_locked_tid_bootstrap.csv")
decision_path <- file.path(output_dir, "satsa_locked_decision.csv")
evaluation_manifest_path <- file.path(output_dir, "satsa_locked_evaluation_manifest.csv")
plot_path <- file.path(output_dir, "satsa_locked_diagnostics.png")
result_paths <- c(metric_path, prediction_path, bootstrap_path, decision_path,
                  evaluation_manifest_path, plot_path)

for (path in c(policy_path, candidate_model_path, standard_model_path, candidate_lock_path)) {
  if (!file.exists(path)) stop("필수 잠금 입력 파일이 없습니다: ", path)
}
if (!requireNamespace("glmnet", quietly = TRUE)) stop("잠금 glmnet 모델을 읽으려면 glmnet 패키지가 필요합니다.")

policy <- read.csv(policy_path, stringsAsFactors = FALSE, check.names = FALSE)
required_policy_columns <- c(
  "policy_version", "dataset", "comparison", "older_age_min",
  "overall_mae_noninferiority_margin", "older_mae_improvement_min",
  "absolute_bias_worsening_margin", "bootstrap_replicates", "bootstrap_seed",
  "confidence_level", "cluster_unit"
)
missing_policy_columns <- setdiff(required_policy_columns, names(policy))
if (length(missing_policy_columns) || nrow(policy) != 1L) {
  stop("SATSA 평가 정책은 필수 열을 가진 한 행이어야 합니다: ",
       paste(missing_policy_columns, collapse = ", "))
}
if (policy$dataset[[1]] != "S-BSST1206" || policy$cluster_unit[[1]] != "TID") {
  stop("SATSA 데이터셋 또는 cluster 단위 정책이 예상과 다릅니다.")
}
numeric_policy <- c(
  policy$older_age_min[[1]], policy$overall_mae_noninferiority_margin[[1]],
  policy$older_mae_improvement_min[[1]], policy$absolute_bias_worsening_margin[[1]],
  policy$bootstrap_replicates[[1]], policy$bootstrap_seed[[1]], policy$confidence_level[[1]]
)
if (any(!is.finite(numeric_policy)) || policy$bootstrap_replicates[[1]] < 1000L ||
    policy$confidence_level[[1]] <= 0 || policy$confidence_level[[1]] >= 1) {
  stop("SATSA 평가 정책의 숫자 값이 올바르지 않습니다.")
}

candidate_model <- readRDS(candidate_model_path)
standard_model <- readRDS(standard_model_path)
candidate_features <- as.character(candidate_model$feature_names)
standard_features <- as.character(standard_model$feature_names)
if (length(candidate_features) != expected_features || anyDuplicated(candidate_features) ||
    !identical(candidate_features, standard_features)) {
  stop("두 잠금 모델은 동일한 순서의 고유한 869개 CpG를 사용해야 합니다.")
}
candidate_lock <- read.csv(candidate_lock_path, stringsAsFactors = FALSE, check.names = FALSE)
candidate_md5 <- unname(tools::md5sum(candidate_model_path))
standard_md5 <- unname(tools::md5sum(standard_model_path))
if (nrow(candidate_lock) != 1L || candidate_lock$model_md5[[1]] != candidate_md5) {
  stop("age-density 후보 모델과 잠금 manifest의 MD5가 일치하지 않습니다.")
}

cat("평가 정책:", policy$policy_version[[1]], "\n")
cat("비교:", policy$comparison[[1]], "\n")
cat("전체 MAE 비열등 허용:", policy$overall_mae_noninferiority_margin[[1]], "년\n")
cat(policy$older_age_min[[1]], "세 이상 MAE 최소 개선:",
    policy$older_mae_improvement_min[[1]], "년\n")
cat("절대 평균 편향 악화 허용:", policy$absolute_bias_worsening_margin[[1]], "년\n")
cat("TID paired cluster bootstrap:", policy$bootstrap_replicates[[1]], "회 / seed",
    policy$bootstrap_seed[[1]], "\n")
cat("후보 모델 MD5:", candidate_md5, "\n")
cat("표준 모델 MD5:", standard_md5, "\n")

if (dry_run) {
  cat("dry-run: SATSA 외부 bundle을 읽거나 성능을 계산하지 않았습니다.\n")
  quit(status = 0)
}
if (!confirm_external) {
  stop("실제 SATSA 평가는 --confirm-external-evaluation을 명시해야 합니다.")
}
if (force && nchar(trimws(rerun_reason)) < 10L) {
  stop("재평가에는 --force와 10자 이상의 --rerun-reason=... 근거가 필요합니다.")
}
existing_results <- result_paths[file.exists(result_paths)]
if (length(existing_results) && !force) {
  stop("SATSA 외부 평가 결과가 이미 있어 재실행을 차단했습니다: ",
       paste(existing_results, collapse = ", "),
       "\n정당한 재평가일 때만 --force --rerun-reason=... 을 사용하세요.")
}
if (!file.exists(bundle_path)) stop("품질 게이트를 통과한 SATSA bundle이 없습니다: ", bundle_path)

bundle <- readRDS(bundle_path)
required_bundle_names <- c(
  "x", "y", "sample_id", "participant_id", "family_cluster_id", "feature_names",
  "source_accession", "candidate_model_md5", "standard_model_md5", "external_evaluation_used"
)
missing_bundle_names <- setdiff(required_bundle_names, names(bundle))
if (length(missing_bundle_names)) stop("SATSA bundle 필수 항목이 없습니다: ", paste(missing_bundle_names, collapse = ", "))
x <- bundle$x
y <- as.numeric(bundle$y)
sample_id <- as.character(bundle$sample_id)
participant_id <- as.character(bundle$participant_id)
family_id <- as.character(bundle$family_cluster_id)
if (!is.matrix(x) || nrow(x) != expected_samples || ncol(x) != expected_features ||
    length(y) != expected_samples || length(sample_id) != expected_samples ||
    length(participant_id) != expected_samples || length(family_id) != expected_samples) {
  stop("SATSA 외부 bundle은 447명 x 869 CpG와 같은 길이의 메타데이터여야 합니다.")
}
if (!identical(colnames(x), candidate_features) ||
    !identical(as.character(bundle$feature_names), candidate_features)) {
  stop("SATSA bundle의 CpG 이름 또는 순서가 잠금 모델과 다릅니다.")
}
if (anyDuplicated(sample_id) || anyDuplicated(participant_id) || any(!nzchar(family_id))) {
  stop("SATSA sample·participant ID는 고유하고 family TID는 비어 있지 않아야 합니다.")
}
if (any(!is.finite(y)) || sum(y >= policy$older_age_min[[1]]) != expected_age80plus) {
  stop("SATSA 연령 또는 80세 이상 표본 수가 잠금 manifest와 다릅니다.")
}
if (any(x < 0 | x > 1, na.rm = TRUE)) stop("SATSA beta 값이 0~1 범위를 벗어났습니다.")
if (bundle$source_accession != policy$dataset[[1]] || isTRUE(bundle$external_evaluation_used)) {
  stop("SATSA source 또는 외부 평가 미사용 표기가 올바르지 않습니다.")
}
if (bundle$candidate_model_md5 != candidate_md5 || bundle$standard_model_md5 != standard_md5) {
  stop("SATSA bundle에 잠근 모델의 MD5가 연결되지 않았습니다.")
}

impute_with_locked_medians <- function(matrix_value, medians, model_name) {
  medians <- as.numeric(medians)
  if (length(medians) != ncol(matrix_value) || any(!is.finite(medians))) {
    stop(model_name, "의 잠금 결측 대치값이 CpG 수와 일치하지 않습니다.")
  }
  result <- matrix_value
  missing_positions <- which(is.na(result), arr.ind = TRUE)
  if (nrow(missing_positions)) result[missing_positions] <- medians[missing_positions[, "col"]]
  if (any(!is.finite(result))) stop(model_name, " 적용 후에도 비유한 beta 값이 남았습니다.")
  result
}

predict_locked <- function(locked_model, model_name) {
  x_imputed <- impute_with_locked_medians(x, locked_model$imputation_medians, model_name)
  as.numeric(predict(locked_model$model, newx = x_imputed, s = locked_model$lambda))
}

standard_prediction <- predict_locked(standard_model, "standard")
candidate_prediction <- predict_locked(candidate_model, "age_density")
if (any(!is.finite(standard_prediction)) || any(!is.finite(candidate_prediction))) {
  stop("잠금 모델 예측에 비유한 값이 있습니다.")
}

metric_row <- function(model_name, subgroup, truth, prediction) {
  residual <- prediction - truth
  calibration <- if (length(unique(prediction)) > 1L) stats::lm(truth ~ prediction) else NULL
  residual_model <- if (length(unique(truth)) > 1L) stats::lm(residual ~ truth) else NULL
  data.frame(
    dataset = policy$dataset[[1]], model = model_name, subgroup = subgroup,
    n = length(truth), MAE = mean(abs(residual)), RMSE = sqrt(mean(residual^2)),
    MedAE = median(abs(residual)),
    R2 = if (length(unique(truth)) > 1L) 1 - sum(residual^2) / sum((truth - mean(truth))^2) else NA_real_,
    bias_mean = mean(residual),
    residual_age_slope = if (is.null(residual_model)) NA_real_ else unname(coef(residual_model)[[2]]),
    calibration_intercept = if (is.null(calibration)) NA_real_ else unname(coef(calibration)[[1]]),
    calibration_slope = if (is.null(calibration)) NA_real_ else unname(coef(calibration)[[2]]),
    stringsAsFactors = FALSE
  )
}

older <- y >= policy$older_age_min[[1]]
metrics <- rbind(
  metric_row("standard", "overall", y, standard_prediction),
  metric_row("age_density", "overall", y, candidate_prediction),
  metric_row("standard", "age80plus", y[older], standard_prediction[older]),
  metric_row("age_density", "age80plus", y[older], candidate_prediction[older])
)

standard_overall <- metrics[metrics$model == "standard" & metrics$subgroup == "overall", , drop = FALSE]
candidate_overall <- metrics[metrics$model == "age_density" & metrics$subgroup == "overall", , drop = FALSE]
standard_older <- metrics[metrics$model == "standard" & metrics$subgroup == "age80plus", , drop = FALSE]
candidate_older <- metrics[metrics$model == "age_density" & metrics$subgroup == "age80plus", , drop = FALSE]

# 80세 이상 표본만 대상으로 가족 TID를 재표집한다. 한 TID에 두 쌍둥이가 있으면 함께 반복된다.
older_indices <- which(older)
older_clusters <- unique(family_id[older_indices])
cluster_rows <- split(older_indices, factor(family_id[older_indices], levels = older_clusters))
set.seed(as.integer(policy$bootstrap_seed[[1]]))
bootstrap_replicates <- as.integer(policy$bootstrap_replicates[[1]])
bootstrap_improvement <- numeric(bootstrap_replicates)
for (iteration in seq_len(bootstrap_replicates)) {
  sampled_clusters <- sample(older_clusters, length(older_clusters), replace = TRUE)
  sampled_rows <- unlist(cluster_rows[sampled_clusters], use.names = FALSE)
  standard_mae <- mean(abs(standard_prediction[sampled_rows] - y[sampled_rows]))
  candidate_mae <- mean(abs(candidate_prediction[sampled_rows] - y[sampled_rows]))
  bootstrap_improvement[[iteration]] <- standard_mae - candidate_mae
}
alpha <- 1 - policy$confidence_level[[1]]
bootstrap_ci <- unname(stats::quantile(
  bootstrap_improvement, probs = c(alpha / 2, 1 - alpha / 2), type = 8
))
bootstrap <- data.frame(
  iteration = seq_len(bootstrap_replicates),
  age80plus_mae_improvement = bootstrap_improvement,
  stringsAsFactors = FALSE
)

overall_mae_difference <- candidate_overall$MAE[[1]] - standard_overall$MAE[[1]]
older_mae_improvement <- standard_older$MAE[[1]] - candidate_older$MAE[[1]]
absolute_bias_worsening <- abs(candidate_overall$bias_mean[[1]]) - abs(standard_overall$bias_mean[[1]])
decision <- data.frame(
  policy_version = policy$policy_version[[1]],
  dataset = policy$dataset[[1]],
  selected_samples = length(y),
  age80plus_samples = sum(older),
  family_clusters = length(unique(family_id)),
  age80plus_family_clusters = length(older_clusters),
  overall_mae_difference_candidate_minus_standard = overall_mae_difference,
  overall_mae_noninferiority_margin = policy$overall_mae_noninferiority_margin[[1]],
  overall_mae_pass = overall_mae_difference <= policy$overall_mae_noninferiority_margin[[1]],
  age80plus_mae_improvement_standard_minus_candidate = older_mae_improvement,
  age80plus_mae_improvement_min = policy$older_mae_improvement_min[[1]],
  age80plus_mae_pass = older_mae_improvement >= policy$older_mae_improvement_min[[1]],
  absolute_bias_worsening_candidate_minus_standard = absolute_bias_worsening,
  absolute_bias_worsening_margin = policy$absolute_bias_worsening_margin[[1]],
  absolute_bias_pass = absolute_bias_worsening <= policy$absolute_bias_worsening_margin[[1]],
  bootstrap_improvement_mean = mean(bootstrap_improvement),
  bootstrap_improvement_ci_lower = bootstrap_ci[[1]],
  bootstrap_improvement_ci_upper = bootstrap_ci[[2]],
  bootstrap_ci_pass = bootstrap_ci[[1]] > 0,
  replace_standard_model = FALSE,
  stringsAsFactors = FALSE
)
decision$replace_standard_model <- with(
  decision, overall_mae_pass & age80plus_mae_pass & absolute_bias_pass & bootstrap_ci_pass
)

predictions <- data.frame(
  sample_id = sample_id, participant_id = participant_id, family_cluster_id = family_id,
  age_true = y, age80plus = older,
  standard_prediction = standard_prediction,
  age_density_prediction = candidate_prediction,
  standard_residual = standard_prediction - y,
  age_density_residual = candidate_prediction - y,
  absolute_error_improvement = abs(standard_prediction - y) - abs(candidate_prediction - y),
  stringsAsFactors = FALSE
)
evaluation_manifest <- data.frame(
  policy_version = policy$policy_version[[1]],
  policy_md5 = unname(tools::md5sum(policy_path)),
  bundle_md5 = unname(tools::md5sum(bundle_path)),
  candidate_model_md5 = candidate_md5,
  standard_model_md5 = standard_md5,
  bootstrap_replicates = bootstrap_replicates,
  bootstrap_seed = policy$bootstrap_seed[[1]],
  rerun = force,
  rerun_reason = if (force) rerun_reason else "first_locked_evaluation",
  evaluated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(metrics, metric_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(predictions, prediction_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bootstrap, bootstrap_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(decision, decision_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(evaluation_manifest, evaluation_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")

grDevices::png(plot_path, width = 1400, height = 900, res = 120)
old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(y, standard_prediction, pch = 16, cex = 0.6, col = rgb(0.15, 0.35, 0.75, 0.35),
     xlab = "실제 나이", ylab = "예측 나이", main = "SATSA 표준 잠금 모델")
abline(0, 1, col = "red", lwd = 2)
plot(y, candidate_prediction, pch = 16, cex = 0.6, col = rgb(0.1, 0.6, 0.35, 0.35),
     xlab = "실제 나이", ylab = "예측 나이", main = "SATSA age-density 후보")
abline(0, 1, col = "red", lwd = 2)
plot(y, standard_prediction - y, pch = 16, cex = 0.6, col = rgb(0.15, 0.35, 0.75, 0.35),
     xlab = "실제 나이", ylab = "잔차", main = "표준 모델 잔차")
abline(h = 0, col = "red", lwd = 2)
plot(y, candidate_prediction - y, pch = 16, cex = 0.6, col = rgb(0.1, 0.6, 0.35, 0.35),
     xlab = "실제 나이", ylab = "잔차", main = "age-density 잔차")
abline(h = 0, col = "red", lwd = 2)
par(old_par)
invisible(dev.off())

message("=== SATSA 잠금 외부 평가 완료 ===")
print(metrics, row.names = FALSE)
print(decision, row.names = FALSE)
message("모델 교체 판정: ", if (decision$replace_standard_model[[1]]) "PASS" else "FAIL")
