#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
project_root <- if (length(file_arg) > 0) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

required_files <- c(
  "README.md",
  "R/randomwalk.R",
  "R/compare_network_properties.R",
  "docs/data.md"
)

missing <- required_files[!file.exists(file.path(project_root, required_files))]
if (length(missing) > 0) {
  stop("Missing required files:\n", paste(missing, collapse = "\n"), call. = FALSE)
}

r_files <- list.files(file.path(project_root, "R"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)

parse_one <- function(path) {
  tryCatch(
    {
      parse(path)
      TRUE
    },
    error = function(err) {
      message("Parse failed: ", path)
      message(conditionMessage(err))
      FALSE
    }
  )
}

ok <- vapply(r_files, parse_one, logical(1))
if (!all(ok)) {
  stop("One or more R files failed to parse.", call. = FALSE)
}

message("Project structure and R syntax checks passed.")
