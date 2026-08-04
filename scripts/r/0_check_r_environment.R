#!/usr/bin/env Rscript

# 원천 데이터를 읽지 않고 현행 분석 파이프라인의 R 실행환경을 점검한다.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("scripts/r/0_check_r_environment.R", winslash = "/", mustWork = TRUE)
}
project_dir <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
manifest_path <- file.path(project_dir, "config", "r_environment.csv")

if (!file.exists(manifest_path)) stop("R 환경 명세가 없습니다: ", manifest_path, call. = FALSE)
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c(
  "component", "component_type", "required_version", "version_policy",
  "scope", "required_by", "reason"
)
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns)) {
  stop("R 환경 명세의 필수 열이 없습니다: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (!nrow(manifest) || anyDuplicated(manifest$component) || any(!nzchar(manifest$component))) {
  stop("R 환경 명세의 component는 비어 있지 않고 고유해야 합니다.", call. = FALSE)
}
if (any(!manifest$component_type %in% c("r", "package"))) {
  stop("component_type은 r 또는 package만 허용합니다.", call. = FALSE)
}
if (any(!manifest$version_policy %in% c("exact", "minimum"))) {
  stop("version_policy는 exact 또는 minimum만 허용합니다.", call. = FALSE)
}

installed_version <- function(component, component_type) {
  if (component_type == "r") return(as.character(getRversion()))
  if (!requireNamespace(component, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(component))
}

version_matches <- function(installed, required, policy) {
  if (is.na(installed)) return(FALSE)
  installed <- numeric_version(installed)
  required <- numeric_version(required)
  if (policy == "exact") installed == required else installed >= required
}

audit <- data.frame(
  component = manifest$component,
  required_version = manifest$required_version,
  version_policy = manifest$version_policy,
  installed_version = mapply(
    installed_version, manifest$component, manifest$component_type,
    USE.NAMES = FALSE
  ),
  stringsAsFactors = FALSE
)
audit$status <- mapply(
  version_matches, audit$installed_version, audit$required_version, audit$version_policy,
  USE.NAMES = FALSE
)
audit$result <- ifelse(audit$status, "PASS", "FAIL")

print(audit[c("component", "required_version", "version_policy", "installed_version", "result")],
      row.names = FALSE)
if (!all(audit$status)) {
  failed <- audit$component[!audit$status]
  stop("R 실행환경 점검 실패: ", paste(failed, collapse = ", "), call. = FALSE)
}

cat("현행 R 실행환경 점검을 통과했습니다.\n")
