#!/usr/bin/env Rscript

# 패키지와 원천 데이터를 불러오지 않고 R 스크립트의 구문만 검사한다.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("scripts/ci/validate_r_syntax.R", winslash = "/", mustWork = TRUE)
}
root_dir <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
r_dir <- file.path(root_dir, "scripts", "r")
r_files <- sort(list.files(r_dir, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE))

if (length(r_files) == 0L) {
  stop("scripts/r 아래에서 검사할 R 파일을 찾지 못했습니다.", call. = FALSE)
}

failures <- character()
for (path in r_files) {
  result <- tryCatch(
    {
      parse(file = path, encoding = "UTF-8")
      NULL
    },
    error = function(error) conditionMessage(error)
  )
  if (!is.null(result)) {
    relative_path <- substring(normalizePath(path, winslash = "/"), nchar(root_dir) + 2L)
    failures <- c(failures, sprintf("%s: %s", relative_path, result))
  }
}

if (length(failures) > 0L) {
  writeLines(c("R 스크립트 구문 검사에 실패했습니다:", paste0("- ", failures)), con = stderr())
  quit(status = 1L)
}

cat(sprintf("R 스크립트 %d개의 구문 검사를 통과했습니다.\n", length(r_files)))
