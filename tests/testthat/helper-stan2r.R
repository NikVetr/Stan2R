skip_if_no_stanc <- function() {
  ok <- tryCatch({ Stan2R:::stanc_path(); TRUE }, error = function(e) FALSE)
  if (!ok) skip("stanc (CmdStan) not available")
}

model_path <- function(name) testthat::test_path("models", paste0(name, ".stan"))

# translate Stan source given as one string
tr <- function(code, mode = "prior_predictive", ...) {
  suppressWarnings(stan2r(stan_program(code), mode, ...))
}

# R lines of a generated script, without blank lines and section rules
code_lines <- function(code) {
  x <- strsplit(code, "\n", fixed = TRUE)[[1]]
  x[nzchar(trimws(x)) & !grepl("^#", trimws(x))]
}

# evaluate a generated script with data `dat` and return `out`
run_code <- function(code, dat, seed = 1) {
  env <- new.env()
  env$dat <- dat
  set.seed(seed)
  eval(parse(text = code), envir = env)
  env$out
}

# a minimal program around transformed-data statements, for expression tests
td_program <- function(data_decls, td_stmts) {
  paste0("data {\n", paste0("  ", data_decls, collapse = "\n"), "\n}\ntransformed data {\n",
         paste0("  ", td_stmts, collapse = "\n"), "\n}\n")
}
