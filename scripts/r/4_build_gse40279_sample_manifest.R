#!/usr/bin/env Rscript

# GSE40279 sample key의 beta matrix 열 이름을 GEO 표본 메타데이터의 연령과 연결한다.

args <- commandArgs(trailingOnly = TRUE)
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."
sample_key_path <- file.path(project_dir, "data", "raw", "GSE40279", "GSE40279_sample_key.txt.gz")
metadata_path <- file.path(project_dir, "outputs", "cohort_metadata_audit", "cohort_sample_metadata.csv")
output_path <- file.path(project_dir, "data", "processed", "GSE40279_sample_manifest.csv")

for (path in c(sample_key_path, metadata_path)) {
  if (!file.exists(path)) stop("필수 입력 파일을 찾을 수 없습니다: ", path)
}

sample_key <- read.delim(
  sample_key_path,
  header = FALSE,
  sep = "\t",
  col.names = c("sample_key_index", "participant_id", "array_id"),
  stringsAsFactors = FALSE
)
sample_key$beta_column_name <- paste0("X", sample_key$participant_id)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- metadata[metadata$series_id == "GSE40279", , drop = FALSE]

# GSE40279의 title은 예: "age 67y 1001"이며 마지막 숫자가 표본 키의 participant_id다.
metadata$participant_id <- sub("^.*\\s([0-9]+)$", "\\1", trimws(metadata$title))
metadata$participant_id[!grepl("^[0-9]+$", metadata$participant_id)] <- NA_character_

if (anyDuplicated(sample_key$participant_id)) stop("sample key의 participant_id가 중복됩니다.")
if (anyDuplicated(metadata$participant_id[!is.na(metadata$participant_id)])) stop("GEO 메타데이터의 participant_id가 중복됩니다.")

manifest <- merge(
  sample_key,
  metadata[, c("participant_id", "sample_id", "age_numeric", "sex", "tissue_detected")],
  by = "participant_id",
  all.x = TRUE,
  sort = FALSE
)
manifest <- manifest[order(manifest$sample_key_index), , drop = FALSE]

if (nrow(manifest) != 656L) stop("예상한 656개 표본 키와 일치하지 않습니다: ", nrow(manifest))
if (anyNA(manifest$sample_id)) stop("GEO 표본 메타데이터와 결합되지 않은 표본이 있습니다.")
if (anyNA(manifest$age_numeric)) stop("연령 정보가 없는 표본이 있습니다.")
if (!all(manifest$tissue_detected == "whole blood")) stop("whole blood가 아닌 표본이 포함되었습니다.")
if (anyDuplicated(manifest$beta_column_name)) stop("beta matrix 열 이름이 중복됩니다.")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(manifest, output_path, row.names = FALSE, na = "")
message("표본 매핑 완료: ", output_path, " (", nrow(manifest), "개)")
