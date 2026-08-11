#!/usr/bin/env Rscript

# GEO의 SOFT 메타데이터만 사용해 후보 코호트의 적합성을 점검한다.
# beta value와 IDAT 원천 파일은 이 스크립트에서 내려받지 않는다.

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
project_arg <- args[!grepl("^--", args)]
# Windows R은 한글이 포함된 경로에 normalizePath()를 적용하면 인코딩이 바뀌어
# file.exists()가 실패할 수 있다. 사용자가 입력한 경로 또는 현재 경로를 그대로 유지한다.
project_dir <- if (length(project_arg)) project_arg[[1]] else "."

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

soft_field <- function(lines, field_name, multiple = FALSE) {
  prefix <- paste0("!", field_name, " = ")
  matches <- lines[startsWith(lines, prefix)]
  values <- sub(paste0("^", prefix), "", matches)
  if (multiple) return(values)
  if (length(values)) values[[1]] else NA_character_
}

to_sample_record <- function(sample_lines, series_id) {
  characteristics <- soft_field(sample_lines, "Sample_characteristics_ch1", multiple = TRUE)
  source_name <- soft_field(sample_lines, "Sample_source_name_ch1")
  title <- soft_field(sample_lines, "Sample_title")
  combined <- paste(c(title, source_name, characteristics), collapse = " | ")

  age_raw <- extract_characteristic(characteristics, "^(age|age \\([^)]*\\)|age at).*:")
  if (is.na(age_raw)) age_raw <- title
  sex <- extract_characteristic(characteristics, "^(sex|gender).*:")
  disease_status <- extract_characteristic(characteristics, "^(disease|disease state|diagnosis|case control).*:")
  smoking_status <- extract_characteristic(characteristics, "^(smoking|smoker|tobacco).*:")

  data.frame(
    series_id = series_id,
    sample_id = soft_field(sample_lines, "Sample_geo_accession"),
    platform_id = soft_field(sample_lines, "Sample_platform_id"),
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

split_soft_samples <- function(lines) {
  start_indices <- which(startsWith(lines, "^SAMPLE = "))
  if (!length(start_indices)) return(list())
  end_indices <- c(start_indices[-1] - 1L, length(lines))
  Map(function(start, end) lines[start:end], start_indices, end_indices)
}

fetch_series_metadata <- function(series_id) {
  message("메타데이터 수집: ", series_id)
  # NCBI GEO의 brief 형식은 표본 특성만 제공하며, beta value·IDAT·데이터 테이블을 포함하지 않는다.
  source_url <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", series_id,
    "&targ=gsm&view=brief&form=text"
  )
  cache_path <- file.path(cache_dir, paste0(series_id, "_gsm_brief.soft"))
  if (!file.exists(cache_path)) {
    utils::download.file(source_url, cache_path, method = "libcurl", mode = "wb", quiet = TRUE)
  }
  lines <- readLines(cache_path, warn = FALSE, encoding = "UTF-8")
  samples <- split_soft_samples(lines)
  if (!length(samples)) stop(series_id, "에서 GSM 메타데이터를 찾지 못했습니다.")
  do.call(rbind, lapply(samples, to_sample_record, series_id = series_id))
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
