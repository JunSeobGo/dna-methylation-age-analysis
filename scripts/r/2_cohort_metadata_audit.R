#!/usr/bin/env Rscript

# GEO의 SOFT 메타데이터만 사용해 후보 코호트의 적합성을 점검한다.
# beta value와 IDAT 원천 파일은 이 스크립트에서 내려받지 않는다.

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
project_arg <- args[!grepl("^--", args)]
project_dir <- normalizePath(if (length(project_arg)) project_arg[[1]] else ".")

registry_path <- file.path(project_dir, "config", "cohort_registry.csv")
output_dir <- file.path(project_dir, "outputs", "cohort_metadata_audit")
cache_dir <- file.path(output_dir, "geo_soft_cache")

if (!file.exists(registry_path)) {
  stop("코호트 registry를 찾을 수 없습니다: ", registry_path)
}

registry <- read.csv(registry_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c(
  "series_id", "classification", "proposed_role", "expected_tissue",
  "expected_platform", "include_in_age_model", "reason", "source_url"
)
missing_columns <- setdiff(required_columns, names(registry))
if (length(missing_columns)) {
  stop("registry에 필요한 열이 없습니다: ", paste(missing_columns, collapse = ", "))
}

if (dry_run) {
  message("--dry-run: registry 형식만 점검합니다. GEO 메타데이터는 내려받지 않습니다.")
  print(registry[, required_columns])
  quit(status = 0)
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop(
    "GEOquery 패키지가 필요합니다. R 콘솔에서 다음을 실행하세요:\n",
    "install.packages('BiocManager')\nBiocManager::install('GEOquery')"
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

extract_characteristic <- function(characteristics, key_pattern) {
  hits <- characteristics[grepl(key_pattern, characteristics, ignore.case = TRUE, perl = TRUE)]
  if (!length(hits)) return(NA_character_)
  trimws(sub("^[^:]+:\\s*", "", hits[[1]]))
}

extract_age <- function(value) {
  if (is.na(value)) return(NA_real_)
  number <- regmatches(value, regexpr("[0-9]+(?:\\.[0-9]+)?", value, perl = TRUE))
  if (!length(number) || !nzchar(number)) return(NA_real_)
  as.numeric(number)
}

detect_tissue <- function(text) {
  if (grepl("whole blood", text, ignore.case = TRUE)) return("whole blood")
  if (grepl("peripheral blood", text, ignore.case = TRUE)) return("peripheral blood")
  if (grepl("blood", text, ignore.case = TRUE)) return("blood (unspecified)")
  NA_character_
}

to_sample_record <- function(gsm, series_id) {
  meta <- GEOquery::Meta(gsm)
  characteristics <- unlist(meta[grep("^characteristics_ch1", names(meta))], use.names = FALSE)
  characteristics <- as.character(characteristics)
  source_name <- as.character(unlist(meta[grep("^source_name_ch1", names(meta))], use.names = FALSE))
  title <- as.character(meta$title %||% NA_character_)
  combined <- paste(c(title, source_name, characteristics), collapse = " | ")

  age_raw <- extract_characteristic(characteristics, "^(age|age \\(years\\)|age at).*:")
  sex <- extract_characteristic(characteristics, "^(sex|gender).*:")
  disease_status <- extract_characteristic(characteristics, "^(disease|disease state|diagnosis|case control).*:")
  smoking_status <- extract_characteristic(characteristics, "^(smoking|smoker|tobacco).*:")

  data.frame(
    series_id = series_id,
    sample_id = as.character(meta$geo_accession %||% NA_character_),
    platform_id = as.character(meta$platform_id %||% NA_character_),
    title = title,
    age_raw = age_raw,
    age_numeric = extract_age(age_raw),
    sex = sex,
    tissue_detected = detect_tissue(combined),
    disease_status = disease_status,
    smoking_status = smoking_status,
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(value, fallback) {
  if (is.null(value) || !length(value)) fallback else value
}

fetch_series_metadata <- function(series_id) {
  message("메타데이터 수집: ", series_id)
  gse <- GEOquery::getGEO(
    series_id,
    GSEMatrix = FALSE,
    getGPL = FALSE,
    destdir = cache_dir
  )
  gsm_list <- GEOquery::GSMList(gse)
  if (!length(gsm_list)) {
    stop(series_id, "에서 GSM 메타데이터를 찾지 못했습니다.")
  }
  do.call(rbind, lapply(gsm_list, to_sample_record, series_id = series_id))
}

sample_tables <- lapply(registry$series_id, function(series_id) {
  tryCatch(
    fetch_series_metadata(series_id),
    error = function(error) {
      warning(series_id, " 수집 실패: ", conditionMessage(error))
      data.frame(
        series_id = series_id,
        sample_id = NA_character_, platform_id = NA_character_, title = NA_character_,
        age_raw = NA_character_, age_numeric = NA_real_, sex = NA_character_,
        tissue_detected = NA_character_, disease_status = NA_character_,
        smoking_status = NA_character_, stringsAsFactors = FALSE
      )
    }
  )
})
samples <- do.call(rbind, sample_tables)

summarize_series <- function(series_id) {
  subset <- samples[samples$series_id == series_id, , drop = FALSE]
  expected <- registry[registry$series_id == series_id, , drop = FALSE]
  n_samples <- sum(!is.na(subset$sample_id))
  age_coverage <- if (n_samples) mean(!is.na(subset$age_numeric)) else 0
  expected_platform <- expected$expected_platform[[1]]
  expected_tissue <- expected$expected_tissue[[1]]
  platform_coverage <- if (n_samples) mean(subset$platform_id == expected_platform, na.rm = TRUE) else 0
  tissue_coverage <- if (n_samples) {
    mean(grepl(expected_tissue, subset$tissue_detected, fixed = TRUE), na.rm = TRUE)
  } else 0
  status <- if (expected$classification[[1]] == "exclude") {
    "제외 유지"
  } else if (n_samples > 0 && age_coverage >= 0.95 && platform_coverage >= 0.95 && tissue_coverage >= 0.95) {
    "적합 후보"
  } else {
    "수동 검토 필요"
  }

  data.frame(
    series_id = series_id,
    classification = expected$classification[[1]],
    proposed_role = expected$proposed_role[[1]],
    include_in_age_model = expected$include_in_age_model[[1]],
    n_samples = n_samples,
    n_with_age = sum(!is.na(subset$age_numeric)),
    age_coverage = round(age_coverage, 3),
    min_age = if (all(is.na(subset$age_numeric))) NA_real_ else min(subset$age_numeric, na.rm = TRUE),
    max_age = if (all(is.na(subset$age_numeric))) NA_real_ else max(subset$age_numeric, na.rm = TRUE),
    expected_platform = expected_platform,
    platform_coverage = round(platform_coverage, 3),
    expected_tissue = expected_tissue,
    tissue_coverage = round(tissue_coverage, 3),
    audit_status = status,
    reason = expected$reason[[1]],
    source_url = expected$source_url[[1]],
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, lapply(registry$series_id, summarize_series))
write.csv(samples, file.path(output_dir, "cohort_sample_metadata.csv"), row.names = FALSE, na = "")
write.csv(summary_table, file.path(output_dir, "cohort_audit_summary.csv"), row.names = FALSE, na = "")

message("점검 완료: ", normalizePath(output_dir))
print(summary_table)
