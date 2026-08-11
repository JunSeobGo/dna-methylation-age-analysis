#!/usr/bin/env Rscript

# 잠긴 age_density 후보의 새 고령층 외부 검증을 위해 SATSA 표현형 표본 목록만 고정한다.
# beta 파일은 내려받거나 읽지 않으며, 반복측정은 개인별 가장 이른 450K 표본 하나만 선택한다.

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
download_metadata <- "--download-metadata" %in% args
force <- "--force" %in% args
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."

accession <- "S-BSST1206"
metadata_url <- "https://www.ebi.ac.uk/biostudies/files/S-BSST1206/SATSA_pheno.txt"
metadata_path <- file.path(project_dir, "data", "raw", accession, "SATSA_pheno.txt")
model_path <- file.path(project_dir, "data", "processed", "locked_age_density_candidate.rds")
model_manifest_path <- file.path(project_dir, "outputs", "age_density_candidate_lock_manifest.csv")
train_path <- file.path(project_dir, "data", "processed", "gse207605_training_bundle.rds")
manifest_path <- file.path(project_dir, "data", "processed", "satsa_external_sample_manifest.csv")
summary_path <- file.path(project_dir, "outputs", "satsa_external_manifest_summary.csv")

if (dry_run) {
  cat("외부 후보:", accession, "\n")
  cat("대상 플랫폼: 450K\n")
  cat("표본 규칙: 개인(IID)별 가장 이른 연령 표본 1개\n")
  cat("불확실성 단위: 쌍둥이 가족(TID) cluster\n")
  cat("beta 파일: 다운로드·읽기 안 함\n")
  cat("표현형 다운로드 URL:", metadata_url, "\n")
  quit(status = 0)
}

if (!file.exists(model_path)) {
  stop("age_density 후보 모델을 먼저 잠가야 합니다: ", model_path)
}
if (!file.exists(model_manifest_path)) stop("후보 모델 잠금 manifest가 없습니다: ", model_manifest_path)
model_manifest <- read.csv(model_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(model_manifest) != 1L || !"model_md5" %in% names(model_manifest) ||
    unname(tools::md5sum(model_path)) != model_manifest$model_md5[[1]]) {
  stop("후보 모델과 잠금 manifest의 MD5가 일치하지 않습니다.")
}
if (!file.exists(train_path)) stop("학습 bundle이 없습니다: ", train_path)
if (!force && file.exists(manifest_path)) {
  stop("SATSA 외부 표본 manifest가 이미 잠겨 있습니다: ", manifest_path,
       "\n표본 규칙을 바꿀 근거가 있을 때만 --force를 사용하세요.")
}

if (download_metadata && !file.exists(metadata_path)) {
  dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
  utils::download.file(metadata_url, metadata_path, method = "libcurl", mode = "wb")
}
if (!file.exists(metadata_path)) {
  stop("SATSA 표현형 파일이 없습니다. --download-metadata를 명시해 공개 표현형만 내려받으세요: ", metadata_path)
}

pheno <- read.delim(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("SID", "TID", "IID", "AGE", "SEX", "ZYG", "SMOKE", "CHIP", "STUDY")
missing <- setdiff(required, names(pheno))
if (length(missing)) stop("SATSA 표현형 필수 열이 없습니다: ", paste(missing, collapse = ", "))
if (nrow(pheno) != 1469L || anyDuplicated(pheno$SID)) stop("SATSA 전체 행 수 또는 SID 고유성이 예상과 다릅니다.")
pheno$AGE <- as.numeric(pheno$AGE)
if (any(!is.finite(pheno$AGE))) stop("SATSA AGE에 숫자가 아닌 값이 있습니다.")

pheno450 <- pheno[tolower(pheno$CHIP) == "450k", , drop = FALSE]
if (nrow(pheno450) != 1094L) stop("SATSA 450K 행 수가 예상한 1,094개와 다릅니다.")
ordering <- order(pheno450$IID, pheno450$AGE, pheno450$SID)
pheno450 <- pheno450[ordering, , drop = FALSE]
selected <- pheno450[!duplicated(pheno450$IID), , drop = FALSE]
if (nrow(selected) != 447L || anyDuplicated(selected$IID) || anyDuplicated(selected$SID)) {
  stop("개인별 1개 표본 선택 결과가 예상한 447명과 다릅니다.")
}

train <- readRDS(train_path)
training_ids <- as.character(train$sample_id)
if (length(intersect(selected$SID, training_ids))) stop("SATSA와 학습 bundle의 표본 ID가 중복됩니다.")

manifest <- data.frame(
  external_sample_id = selected$SID,
  participant_id = selected$IID,
  family_cluster_id = selected$TID,
  chronological_age = selected$AGE,
  sex = selected$SEX,
  zygosity = selected$ZYG,
  smoking_status = selected$SMOKE,
  platform = selected$CHIP,
  study_wave = selected$STUDY,
  selection_rule = "earliest_450k_sample_per_IID",
  stringsAsFactors = FALSE
)
summary <- data.frame(
  accession = accession,
  source_rows = nrow(pheno), source_450k_rows = nrow(pheno450),
  selected_samples = nrow(manifest), unique_participants = length(unique(manifest$participant_id)),
  family_clusters = length(unique(manifest$family_cluster_id)),
  age_min = min(manifest$chronological_age), age_max = max(manifest$chronological_age),
  age80plus = sum(manifest$chronological_age >= 80), age90plus = sum(manifest$chronological_age >= 90),
  women = sum(manifest$sex == "Woman"), men = sum(manifest$sex == "Man"),
  training_id_overlap = length(intersect(manifest$external_sample_id, training_ids)),
  beta_data_used = FALSE, model_selection_used = FALSE,
  candidate_model_md5 = model_manifest$model_md5[[1]],
  phenotype_md5 = unname(tools::md5sum(metadata_path)),
  manifest_locked_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)

dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
write.csv(manifest, manifest_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(summary, summary_path, row.names = FALSE, fileEncoding = "UTF-8")

message("=== SATSA 외부 검증 표본 manifest 잠금 완료 ===")
print(summary, row.names = FALSE)
message("beta 파일은 내려받거나 읽지 않았습니다.")
