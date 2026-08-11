#!/usr/bin/env Rscript

# GSE40279의 표본 키와 beta value를 명시적 요청에서만 수집한다.
# 기본 실행은 다운로드 계획만 보여 주므로 대용량 파일이 자동으로 내려받아지지 않는다.

args <- commandArgs(trailingOnly = TRUE)
download_sample_key <- "--download-sample-key" %in% args
download_beta <- "--download-beta" %in% args
force_download <- "--force" %in% args
project_args <- args[!grepl("^--", args)]
project_dir <- if (length(project_args)) project_args[[1]] else "."
manifest_path <- file.path(project_dir, "config", "dataset_download_manifest.csv")

if (!file.exists(manifest_path)) {
  stop("다운로드 매니페스트를 찾을 수 없습니다: ", manifest_path)
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c(
  "dataset_id", "resource_id", "resource_type", "remote_url", "local_path",
  "expected_size_mb", "required_for_stage", "description"
)
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns)) {
  stop("매니페스트에 필요한 열이 없습니다: ", paste(missing_columns, collapse = ", "))
}

gse40279 <- manifest[manifest$dataset_id == "GSE40279", , drop = FALSE]
if (nrow(gse40279) != 2L) stop("GSE40279 리소스 정의가 올바르지 않습니다.")

print(gse40279[, c("resource_id", "expected_size_mb", "local_path", "description")])

if (!download_sample_key && !download_beta) {
  message(
    "다운로드 계획만 표시했습니다. 표본 키만 받으려면 --download-sample-key, ",
    "약 1.1GB beta matrix까지 받으려면 --download-beta를 사용하세요."
  )
  quit(status = 0)
}

is_gzip_file <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  identical(readBin(connection, what = "raw", n = 2L), as.raw(c(0x1f, 0x8b)))
}

download_resource <- function(resource) {
  destination <- file.path(project_dir, resource$local_path[[1]])
  partial_destination <- paste0(destination, ".part")
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(destination) && !force_download) {
    if (is_gzip_file(destination)) {
      message("기존 파일 사용: ", destination)
      return(destination)
    }
    stop("기존 파일의 gzip 형식이 올바르지 않습니다. --force로 다시 받으세요: ", destination)
  }

  if (force_download && file.exists(partial_destination)) unlink(partial_destination)
  partial_size_mb <- if (file.exists(partial_destination)) round(file.info(partial_destination)$size / 1024^2, 1) else 0
  message("다운로드 시작: ", resource$resource_id[[1]], " (약 ", resource$expected_size_mb[[1]], "MB, 현재 ", partial_size_mb, "MB)")

  # Windows curl.exe는 HTTP Range 요청으로 .part 파일부터 재개하며, 네트워크 끊김도 재시도한다.
  curl_path <- Sys.which("curl.exe")
  if (!nzchar(curl_path)) stop("재개 다운로드에 필요한 curl.exe를 찾을 수 없습니다.")
  curl_args <- c(
    "--fail", "--location", "--continue-at", "-", "--retry", "5", "--retry-delay", "5",
    "--output", partial_destination, resource$remote_url[[1]]
  )
  exit_status <- system2(curl_path, curl_args)
  if (!identical(exit_status, 0L)) stop("다운로드가 완료되지 않았습니다. 같은 명령을 다시 실행하면 .part 파일부터 재개합니다.")

  if (!file.exists(partial_destination) || !is_gzip_file(partial_destination)) {
    stop("다운로드 파일의 gzip 형식을 확인하지 못했습니다: ", partial_destination)
  }
  if (file.exists(destination)) unlink(destination)
  if (!file.rename(partial_destination, destination)) {
    stop("임시 파일 이름 변경에 실패했습니다: ", partial_destination)
  }
  message("다운로드 완료: ", destination, " (", round(file.info(destination)$size / 1024^2, 1), "MB)")
  destination
}

if (download_beta) download_sample_key <- TRUE
if (download_sample_key) {
  sample_key <- gse40279[gse40279$resource_id == "sample_key", , drop = FALSE]
  download_resource(sample_key)
}
if (download_beta) {
  beta_matrix <- gse40279[gse40279$resource_id == "average_beta", , drop = FALSE]
  download_resource(beta_matrix)
}
