#!/usr/bin/env Rscript

# 최종 잠금 외부 검증 세트 GSE87571의 beta 행렬에서 869개 공통 CpG를 학습 bundle과
# 동일한 feature 순서로 추출해 외부 검증용 RDS를 준비한다.
#
# 매우 중요
# - GSE87571은 학습·특징 선택·하이퍼파라미터 선택에 절대 사용하지 않는다.
# - 이 단계에서는 모델 성능(나이 기반 지표)을 계산하지 않는다. 행렬 생성과 표본·연령
#   매핑 검증까지만 수행한다.
#
# 데이터 구조
# - GSE87571_matrix1of2.txt.gz, GSE87571_matrix2of2.txt.gz (SWAN 정규화 Average Beta)
# - 표본은 두 파일에 나뉘어 있고(표본 단위 분할), 각 표본은 두 열로 저장된다.
#   Xn   = Average Beta
#   Xn.1 = detection p-value
# - 열 라벨 Xn은 SOFT의 sample_title(예: "X1 genomic DNA from whole blood")과 연결되고,
#   SOFT에서 GSM 접근번호와 age를 얻는다.

args <- commandArgs(trailingOnly = TRUE)
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

raw_dir <- file.path(project_dir, "data", "raw", "GSE87571")
matrix_paths <- file.path(raw_dir, c("GSE87571_matrix1of2.txt.gz", "GSE87571_matrix2of2.txt.gz"))
soft_path <- file.path(project_dir, "outputs", "cohort_metadata_audit", "geo_soft_cache", "GSE87571_gsm_brief.soft")
training_bundle_path <- file.path(project_dir, "data", "processed", "gse207605_training_bundle.rds")
probe_path <- file.path(project_dir, "data", "processed", "gse207605_primary_common_probes.txt")
bundle_path <- file.path(project_dir, "data", "processed", "gse87571_external_validation_bundle.rds")
summary_path <- file.path(project_dir, "outputs", "gse87571_external_validation_summary.csv")

for (p in c(matrix_paths, soft_path, training_bundle_path, probe_path)) {
  if (!file.exists(p)) stop("필요한 입력 파일이 없습니다: ", p)
}

# --- 학습 bundle의 feature 순서 확보 -----------------------------------------

training_bundle <- readRDS(training_bundle_path)
feature_names <- training_bundle$feature_names
training_ids <- training_bundle$sample_id
if (length(feature_names) != 869L) stop("학습 feature 수가 869가 아닙니다: ", length(feature_names))

probe_set <- feature_names  # 학습과 동일 순서로 정렬하기 위한 기준

# --- SOFT 메타데이터 파싱 (title 토큰 -> GSM, age) ---------------------------

parse_soft <- function(path) {
  lines <- readLines(path, warn = FALSE)
  gsm <- token <- age <- disease <- character(0)
  cur_gsm <- cur_token <- cur_age <- cur_disease <- NA_character_
  flush <- function() {
    if (!is.na(cur_gsm)) {
      gsm[[length(gsm) + 1L]] <<- cur_gsm
      token[[length(token) + 1L]] <<- cur_token
      age[[length(age) + 1L]] <<- cur_age
      disease[[length(disease) + 1L]] <<- cur_disease
    }
  }
  for (ln in lines) {
    if (startsWith(ln, "^SAMPLE = ")) {
      flush()
      cur_gsm <- sub("^\\^SAMPLE = ", "", ln)
      cur_token <- cur_age <- cur_disease <- NA_character_
    } else if (startsWith(ln, "!Sample_title = ")) {
      title <- sub("^!Sample_title = ", "", ln)
      cur_token <- strsplit(title, "\\s+")[[1]][1]
    } else if (grepl("!Sample_characteristics_ch1 = age:", ln)) {
      cur_age <- trimws(sub(".*age:\\s*", "", ln))
    } else if (grepl("!Sample_characteristics_ch1 = disease state:", ln)) {
      cur_disease <- trimws(sub(".*disease state:\\s*", "", ln))
    }
  }
  flush()
  data.frame(
    gsm = unlist(gsm), token = unlist(token),
    age = suppressWarnings(as.numeric(unlist(age))),
    disease = unlist(disease), stringsAsFactors = FALSE
  )
}

meta <- parse_soft(soft_path)
if (anyDuplicated(meta$token)) stop("SOFT title 토큰에 중복이 있습니다.")
if (anyDuplicated(meta$gsm)) stop("SOFT GSM에 중복이 있습니다.")
message("SOFT 메타데이터 표본 ", nrow(meta), "개 (연령 유효 ", sum(!is.na(meta$age)), "개)")

# --- beta 행렬 스트리밍 추출 (869 probe만) -----------------------------------

# 큰 gz 파일을 청크로 읽어 첫 열(probe ID)이 목표 집합에 속하는 행만 남긴다(저메모리).
filter_matrix <- function(path, probes) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L)
  col_names <- strsplit(header, "\t", fixed = TRUE)[[1]]
  probe_lookup <- setNames(rep(TRUE, length(probes)), probes)
  kept <- character(0)
  repeat {
    chunk <- readLines(con, n = 50000L)
    if (!length(chunk)) break
    first_token <- sub("\t.*$", "", chunk)
    hit <- !is.na(probe_lookup[first_token])
    if (any(hit)) kept <- c(kept, chunk[hit])
  }
  # 남긴 행을 행렬로 변환한다.
  parts <- strsplit(kept, "\t", fixed = TRUE)
  ids <- vapply(parts, `[`, character(1), 1L)
  values <- t(vapply(parts, function(p) suppressWarnings(as.numeric(p[-1])), numeric(length(col_names) - 1L)))
  rownames(values) <- ids
  colnames(values) <- col_names[-1]
  values
}

message("matrix1of2 추출 중...")
m1 <- filter_matrix(matrix_paths[1], probe_set)
message("  ", nrow(m1), "개 probe × ", ncol(m1), "개 열")
message("matrix2of2 추출 중...")
m2 <- filter_matrix(matrix_paths[2], probe_set)
message("  ", nrow(m2), "개 probe × ", ncol(m2), "개 열")

# --- beta 열과 detection p-value 열 구분 -------------------------------------

# Xn(기본)과 Xn.1(중복) 열을 나눈다. 어느 쪽이 beta인지 exact-zero 비율로 판정한다.
# detection p-value 열은 정확히 0인 값이 매우 많고, beta 열은 그렇지 않다.
split_beta_columns <- function(mat) {
  base_cols <- grep("^X[0-9]+$", colnames(mat), value = TRUE)
  dot_cols <- grep("^X[0-9]+\\.1$", colnames(mat), value = TRUE)
  if (length(base_cols) != length(dot_cols)) {
    stop("beta/pval 열 쌍이 맞지 않습니다: base=", length(base_cols), ", dot=", length(dot_cols))
  }
  zero_frac <- function(cols) {
    sub <- mat[, cols, drop = FALSE]
    mean(sub == 0, na.rm = TRUE)
  }
  base_zero <- zero_frac(base_cols)
  dot_zero <- zero_frac(dot_cols)
  message(sprintf("  exact-zero 비율: base=%.4f, .1=%.4f", base_zero, dot_zero))
  # exact-zero가 더 적은 쪽을 beta로 선택한다.
  if (base_zero <= dot_zero) {
    beta_cols <- base_cols
    tokens <- beta_cols
  } else {
    beta_cols <- dot_cols
    tokens <- sub("\\.1$", "", beta_cols)
  }
  if (!identical(sort(beta_cols), sort(base_cols))) {
    warning("beta 열이 기본(Xn) 열이 아닌 것으로 판정되었습니다. 데이터 형식을 재확인하세요.")
  }
  beta <- mat[, beta_cols, drop = FALSE]
  colnames(beta) <- tokens
  beta
}

message("matrix1of2 beta 열 판정")
b1 <- split_beta_columns(m1)
message("matrix2of2 beta 열 판정")
b2 <- split_beta_columns(m2)

# --- 두 파일 결합 및 정렬 ----------------------------------------------------

# 두 파일은 표본 단위로 분할되어 있으므로 열(표본)을 합친다. probe 행은 학습 순서로 정렬한다.
if (!all(probe_set %in% rownames(b1)) || !all(probe_set %in% rownames(b2))) {
  missing1 <- setdiff(probe_set, rownames(b1))
  missing2 <- setdiff(probe_set, rownames(b2))
  stop("일부 공통 CpG가 matrix 파일에 없습니다. 누락 예: ",
       paste(head(c(missing1, missing2), 3), collapse = ", "))
}
b1 <- b1[probe_set, , drop = FALSE]
b2 <- b2[probe_set, , drop = FALSE]

common_tokens <- intersect(colnames(b1), colnames(b2))
if (length(common_tokens)) stop("두 matrix 파일에 겹치는 표본 토큰이 있습니다: ", length(common_tokens))

beta_all <- cbind(b1, b2)  # probe × 표본
tokens <- colnames(beta_all)

# --- 토큰 -> GSM/age 매핑 ----------------------------------------------------

token_to_gsm <- setNames(meta$gsm, meta$token)
token_to_age <- setNames(meta$age, meta$token)
if (!all(tokens %in% meta$token)) {
  stop("matrix 표본 토큰 중 SOFT에 없는 것이 있습니다: ",
       paste(head(setdiff(tokens, meta$token), 3), collapse = ", "))
}
gsm_ids <- unname(token_to_gsm[tokens])
ages <- unname(token_to_age[tokens])

# --- 표본 × feature 행렬 구성 ------------------------------------------------

x <- t(beta_all)            # 표본 × probe
rownames(x) <- gsm_ids
colnames(x) <- probe_set    # 학습 feature 순서와 동일

# --- 검증 (성능 계산 없음) ---------------------------------------------------

stopifnot(
  "feature 순서가 학습 bundle과 다릅니다" = identical(colnames(x), feature_names),
  "표본 수와 age 길이가 다릅니다" = nrow(x) == length(ages),
  "표본 수와 GSM 길이가 다릅니다" = nrow(x) == length(gsm_ids)
)
if (anyDuplicated(gsm_ids)) stop("외부 검증 표본 GSM에 중복이 있습니다.")

# 학습 데이터와 GSM 중복이 없어야 한다(누수 방지).
leaked <- intersect(gsm_ids, training_ids)
if (length(leaked)) stop("학습 데이터와 GSM이 겹칩니다(누수): ", paste(head(leaked, 5), collapse = ", "))

# beta 범위 검증(NA 허용).
finite_x <- x[is.finite(x)]
if (length(finite_x) && (min(finite_x) < 0 || max(finite_x) > 1)) {
  stop("0~1 범위를 벗어난 beta 값이 있습니다: min=", min(finite_x), ", max=", max(finite_x))
}

na_total <- sum(is.na(x))
n_valid_age <- sum(!is.na(ages))

# --- bundle 저장 -------------------------------------------------------------

bundle <- list(
  x = x,
  y = ages,
  group = rep("GSE87571", nrow(x)),
  sample_id = gsm_ids,
  feature_names = probe_set
)
dir.create(dirname(bundle_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(bundle, bundle_path)

# --- 요약 저장 (나이 분포는 데이터 기술일 뿐 모델 성능이 아님) ---------------

summary_table <- data.frame(
  dataset = "GSE87571",
  samples = nrow(x),
  features = ncol(x),
  valid_age = n_valid_age,
  age_min = if (n_valid_age) min(ages, na.rm = TRUE) else NA_real_,
  age_max = if (n_valid_age) max(ages, na.rm = TRUE) else NA_real_,
  na_cells = na_total,
  na_ratio = na_total / length(x),
  gsm_overlap_with_training = length(leaked),
  feature_order_matches_training = identical(colnames(x), feature_names),
  stringsAsFactors = FALSE
)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, summary_path, row.names = FALSE, fileEncoding = "UTF-8")

message("\n외부 검증 bundle 저장 완료: ", bundle_path)
message(sprintf("X: %d개 표본 × %d개 CpG, 연령 유효 %d개, 결측 %d셀", nrow(x), ncol(x), n_valid_age, na_total))
print(summary_table, row.names = FALSE)
message("주의: 이 단계에서는 모델 성능을 계산하지 않았습니다.")
