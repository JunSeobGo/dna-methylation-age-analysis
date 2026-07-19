#!/usr/bin/env Rscript

# 1차 학습 bundle을 사용해 코호트 그룹 기반 중첩 교차검증으로 Elastic Net 벤치마크를 학습·평가한다.
#
# 데이터 누수 방지 설계
# - outer loop: leave-one-cohort-out. outer test 코호트는 inner tuning에 전혀 쓰지 않는다.
# - inner loop: outer 학습 코호트 안에서 다시 leave-one-cohort-out 그룹 CV.
# - 결측 대치(중앙값)와 표준화는 각 학습 fold 안에서만 계산하고, 검증·테스트 fold에는
#   학습 fold에서 계산한 값만 적용한다. 전체 데이터 기준 대치·표준화를 하지 않는다.
# - glmnet의 내부 표준화(standardize=TRUE)는 fit에 넘긴 학습 fold 데이터만으로 계산되므로
#   누수가 없다.
#
# GSE55763(2,639명)이 표본 수로 모델을 지배하지 않도록 두 가중치 방식을 비교한다.
# - standard: 모든 표본 가중치 1
# - balanced: 코호트별 총 가중치가 같도록 표본 가중치 = 1 / (해당 코호트 표본 수), 합이 n이 되게 정규화
#
# 주의: 869개 공통 CpG는 GSE207605 원 논문에서 연령 연관성으로 사전 선택된 특징이므로 이 모델은
# 벤치마크로 해석한다. 최종 일반화 성능은 GSE87571 외부 검증(별도 단계)에서 판단한다.

suppressPackageStartupMessages(library(glmnet))
set.seed(20260719)

args <- commandArgs(trailingOnly = TRUE)
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

bundle_path <- file.path(project_dir, "data", "processed", "gse207605_training_bundle.rds")
out_dir <- file.path(project_dir, "outputs")
if (!file.exists(bundle_path)) {
  stop("학습 bundle을 찾을 수 없습니다: ", bundle_path, "\n먼저 7_build_gse207605_training_bundle.R을 실행하세요.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(bundle_path)
x <- bundle$x
y <- bundle$y
group <- bundle$group
sample_id <- bundle$sample_id

stopifnot(
  "x가 행렬이 아닙니다" = is.matrix(x),
  "행 수 불일치" = nrow(x) == length(y),
  "group 길이 불일치" = length(group) == length(y),
  "sample_id 길이 불일치" = length(sample_id) == length(y)
)

alphas <- c(0.05, 0.1, 0.25, 0.5, 0.75, 1.0)
cohorts <- sort(unique(group))
message("코호트 ", length(cohorts), "개로 leave-one-cohort-out 중첩 교차검증을 수행합니다.")

# --- 보조 함수 ---------------------------------------------------------------

# 학습 fold의 열별 중앙값(결측 제외). 열 전체가 NA이면 0.5로 대체(방어적, 실제로는 발생하지 않음).
column_medians <- function(mat) {
  med <- apply(mat, 2, function(col) median(col, na.rm = TRUE))
  med[is.na(med)] <- 0.5
  med
}

# 학습 fold에서 계산한 중앙값으로만 결측 대치.
impute_with <- function(mat, med) {
  na_pos <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(na_pos)) mat[na_pos] <- med[na_pos[, "col"]]
  mat
}

# 가중치 계산. balanced는 코호트별 총 가중치가 같아지도록 하고 합이 n이 되게 정규화한다.
make_weights <- function(g, scheme) {
  if (scheme == "standard") return(rep(1, length(g)))
  sizes <- table(g)
  w <- as.numeric(1 / sizes[as.character(g)])
  w * length(g) / sum(w)
}

# --- 중첩 교차검증 (한 가중치 방식) -----------------------------------------

run_scheme <- function(scheme) {
  message("== 가중치 방식: ", scheme, " ==")
  oof_pred <- rep(NA_real_, length(y))
  fold_rows <- list()

  for (o in cohorts) {
    test_idx <- which(group == o)
    train_idx <- which(group != o)
    xtr <- x[train_idx, , drop = FALSE]
    ytr <- y[train_idx]
    gtr <- group[train_idx]
    xte <- x[test_idx, , drop = FALSE]
    train_cohorts <- sort(unique(gtr))

    # outer 학습 fold의 중앙값으로 학습·테스트를 대치(테스트는 학습 값만 사용).
    med_outer <- column_medians(xtr)
    xtr_imp <- impute_with(xtr, med_outer)
    xte_imp <- impute_with(xte, med_outer)
    w_outer <- make_weights(gtr, scheme)

    best <- list(mae = Inf, alpha = NA_real_, lambda = NA_real_)

    for (a in alphas) {
      # 후보 lambda 격자는 outer 학습 fit에서 유도해 모든 inner fold에 동일하게 적용한다.
      lambda_grid <- glmnet(xtr_imp, ytr, alpha = a, weights = w_outer, standardize = TRUE)$lambda
      inner_pred <- matrix(NA_real_, nrow = length(ytr), ncol = length(lambda_grid))

      for (ic in train_cohorts) {
        iv <- which(gtr == ic)
        it <- which(gtr != ic)
        med_in <- column_medians(xtr[it, , drop = FALSE])
        xit <- impute_with(xtr[it, , drop = FALSE], med_in)
        xiv <- impute_with(xtr[iv, , drop = FALSE], med_in)
        w_in <- make_weights(gtr[it], scheme)
        fit <- glmnet(xit, ytr[it], alpha = a, lambda = lambda_grid, weights = w_in, standardize = TRUE)
        pr <- predict(fit, newx = xiv)
        inner_pred[iv, seq_len(ncol(pr))] <- pr
      }

      # inner CV MAE(표본 단위 평균)로 lambda를 고른다. 모든 inner val 표본이 채워졌을 때만 유효.
      valid <- colSums(is.na(inner_pred)) == 0
      if (!any(valid)) next
      mae_lambda <- rep(Inf, length(lambda_grid))
      mae_lambda[valid] <- colMeans(abs(inner_pred[, valid, drop = FALSE] - ytr))
      k <- which.min(mae_lambda)
      if (mae_lambda[k] < best$mae) {
        best <- list(mae = mae_lambda[k], alpha = a, lambda = lambda_grid[k])
      }
    }

    # 최종 refit: outer 학습 전체에 best alpha, best lambda로 학습해 held-out 코호트를 예측한다.
    fit <- glmnet(xtr_imp, ytr, alpha = best$alpha, weights = w_outer, standardize = TRUE)
    pr <- predict(fit, newx = xte_imp, s = best$lambda)
    oof_pred[test_idx] <- as.numeric(pr)

    coefs <- as.numeric(coef(fit, s = best$lambda))[-1]
    nonzero <- sum(coefs != 0)
    fold_rows[[o]] <- data.frame(
      scheme = scheme, outer_cohort = o, n_test = length(test_idx),
      alpha = best$alpha, lambda = best$lambda, inner_cv_mae = best$mae,
      nonzero = nonzero, stringsAsFactors = FALSE
    )
    message(sprintf(
      "  [%s] test=%s n=%d alpha=%.2f lambda=%.4f nonzero=%d",
      scheme, o, length(test_idx), best$alpha, best$lambda, nonzero
    ))
  }

  list(pred = oof_pred, params = do.call(rbind, fold_rows))
}

# --- 지표 계산 ---------------------------------------------------------------

metrics_overall <- function(scheme, pred) {
  err <- pred - y
  ss_res <- sum(err^2)
  ss_tot <- sum((y - mean(y))^2)
  data.frame(
    scheme = scheme,
    n = length(y),
    MAE = mean(abs(err)),
    RMSE = sqrt(mean(err^2)),
    MedAE = median(abs(err)),
    R2 = 1 - ss_res / ss_tot,
    bias_mean = mean(err),
    resid_age_slope = as.numeric(coef(lm(err ~ y))[2]),
    stringsAsFactors = FALSE
  )
}

metrics_by_cohort <- function(scheme, pred) {
  do.call(rbind, lapply(cohorts, function(c) {
    idx <- which(group == c)
    err <- pred[idx] - y[idx]
    data.frame(
      scheme = scheme, cohort = c, n = length(idx),
      age_min = min(y[idx]), age_max = max(y[idx]),
      MAE = mean(abs(err)), RMSE = sqrt(mean(err^2)),
      MedAE = median(abs(err)), bias_mean = mean(err),
      stringsAsFactors = FALSE
    )
  }))
}

age_breaks <- c(-Inf, 20, 40, 60, 80, Inf)
age_labels <- c("<20", "20-39", "40-59", "60-79", ">=80")
metrics_by_agebin <- function(scheme, pred) {
  bin <- cut(y, breaks = age_breaks, labels = age_labels, right = FALSE)
  do.call(rbind, lapply(levels(bin), function(b) {
    idx <- which(bin == b)
    if (!length(idx)) return(NULL)
    err <- pred[idx] - y[idx]
    data.frame(
      scheme = scheme, age_bin = b, n = length(idx),
      MAE = mean(abs(err)), RMSE = sqrt(mean(err^2)),
      MedAE = median(abs(err)), bias_mean = mean(err),
      stringsAsFactors = FALSE
    )
  }))
}

# --- 실행 --------------------------------------------------------------------

schemes <- c("standard", "balanced")
results <- lapply(schemes, run_scheme)
names(results) <- schemes

overall <- do.call(rbind, lapply(schemes, function(s) metrics_overall(s, results[[s]]$pred)))
by_cohort <- do.call(rbind, lapply(schemes, function(s) metrics_by_cohort(s, results[[s]]$pred)))
by_agebin <- do.call(rbind, lapply(schemes, function(s) metrics_by_agebin(s, results[[s]]$pred)))
fold_params <- do.call(rbind, lapply(schemes, function(s) results[[s]]$params))

oof <- data.frame(
  sample_id = sample_id,
  cohort = group,
  age_true = y,
  pred_standard = results[["standard"]]$pred,
  pred_balanced = results[["balanced"]]$pred,
  stringsAsFactors = FALSE
)

write.csv(overall, file.path(out_dir, "gse207605_cv_metrics_overall.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(by_cohort, file.path(out_dir, "gse207605_cv_metrics_by_cohort.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(by_agebin, file.path(out_dir, "gse207605_cv_metrics_by_agebin.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fold_params, file.path(out_dir, "gse207605_cv_fold_params.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(oof, file.path(out_dir, "gse207605_cv_oof_predictions.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# --- 진단 그래프 -------------------------------------------------------------

png(file.path(out_dir, "gse207605_cv_diagnostics.png"), width = 1200, height = 1000, res = 120)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (s in schemes) {
  pred <- results[[s]]$pred
  plot(y, pred, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.25),
       xlab = "실제 나이", ylab = "예측 나이",
       main = sprintf("%s: 예측 vs 실제 (MAE %.2f)", s, mean(abs(pred - y))))
  abline(0, 1, col = "red", lwd = 2)
  plot(y, pred - y, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.25),
       xlab = "실제 나이", ylab = "잔차 (예측 - 실제)",
       main = sprintf("%s: 잔차 vs 실제", s))
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(y, pred - y), col = "blue", lwd = 2)
}
par(op)
invisible(dev.off())

# --- 콘솔 요약 ---------------------------------------------------------------

message("\n=== 전체 OOF 지표 ===")
print(overall, row.names = FALSE)
message("\n=== 코호트별 지표 ===")
print(by_cohort, row.names = FALSE)
message("\n=== 연령 구간별 지표 ===")
print(by_agebin, row.names = FALSE)
message("\n=== fold별 선택 파라미터 ===")
print(fold_params, row.names = FALSE)
message("\n산출물을 outputs/에 저장했습니다(Git 제외).")
