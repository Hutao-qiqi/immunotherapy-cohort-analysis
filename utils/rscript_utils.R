get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[grepl("^--file=", cmd)]
  if (length(file_arg) == 0) {
    return(NA_character_)
  }
  sub("^--file=", "", file_arg[1])
}

get_script_dir <- function() {
  script_path <- get_script_path()
  if (is.na(script_path) || !nzchar(script_path)) {
    return(getwd())
  }
  normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
}

# Assumes the repository/workspace root is 4 levels above:
# iscience/codes/github/<module>/<script>.R
get_project_root <- function(script_dir = get_script_dir(), levels_up = 4) {
  stopifnot(is.numeric(levels_up), levels_up >= 0)
  parts <- rep("..", levels_up)
  normalizePath(do.call(file.path, c(list(script_dir), as.list(parts))), winslash = "/", mustWork = FALSE)
}

set_working_dir_to_script <- function() {
  wd <- get_script_dir()
  setwd(wd)
  invisible(wd)
}
