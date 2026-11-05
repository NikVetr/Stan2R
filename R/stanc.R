# IO boundary to the Stan compiler. Everything else in Stan2R works on text that
# stanc printed, so a browser port only has to replace this file with calls to
# stanc.js (`stanc(name, code, ["debug-typed-ast"])`).

stanc_path <- function() {
  opt <- getOption("stan2r.stanc", Sys.getenv("STAN2R_STANC", ""))
  if (nzchar(opt)) {
    if (!file.exists(opt)) stop("Stan2R: stanc not found at '", opt, "' (option stan2r.stanc / STAN2R_STANC)", call. = FALSE)
    return(normalizePath(opt))
  }
  candidates <- character(0)
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    cmdstan <- tryCatch(suppressMessages(cmdstanr::cmdstan_path()), error = function(e) "")
    if (nzchar(cmdstan)) candidates <- c(candidates, file.path(cmdstan, "bin", "stanc"))
  }
  installs <- sort(Sys.glob(file.path(path.expand("~"), ".cmdstan", "cmdstan-*", "bin", "stanc")), decreasing = TRUE)
  candidates <- c(candidates, installs)
  if (.Platform$OS.type == "windows") candidates <- c(paste0(candidates, ".exe"), candidates)
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop("Stan2R: could not find the stanc compiler. Install CmdStan (cmdstanr::install_cmdstan()) ",
         "or set options(stan2r.stanc = '/path/to/cmdstan/bin/stanc').", call. = FALSE)
  }
  found[1]
}

stanc_version_string <- function() {
  path <- tryCatch(stanc_path(), error = function(e) NA_character_)
  if (is.na(path)) return("unknown")
  v <- tryCatch(system2(path, "--version", stdout = TRUE, stderr = TRUE), error = function(e) "unknown")
  sub("^stanc3 v?([^ ]+).*$", "\\1", v[1])
}

# Run stanc on Stan source text; returns list(stdout, warnings). A compile error
# stops with stanc's own message, which carries the line and column.
run_stanc <- function(stan_code, args, model_name = "model") {
  dir <- tempfile("stan2r_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  model_file <- file.path(dir, paste0(model_name, ".stan"))
  writeLines(stan_code, model_file, useBytes = TRUE)
  out_file <- file.path(dir, "stdout.txt")
  err_file <- file.path(dir, "stderr.txt")
  status <- system2(stanc_path(),
                    c(args, "--color=never", paste0("--filename-in-msg=", model_name, ".stan"),
                      "-o", shQuote(file.path(dir, "out.hpp")), shQuote(model_file)),
                    stdout = out_file, stderr = err_file)
  stdout <- if (file.exists(out_file)) readLines(out_file, warn = FALSE) else character(0)
  stderr <- if (file.exists(err_file)) readLines(err_file, warn = FALSE) else character(0)
  output <- file.path(dir, "out.hpp")
  output <- if (file.exists(output)) readLines(output, warn = FALSE) else character(0)
  if (!identical(as.integer(status), 0L)) {
    msg <- paste(c(stderr, stdout), collapse = "\n")
    stop("Stan2R: stanc rejected the model:\n", msg, call. = FALSE)
  }
  list(stdout = stdout, output = output, warnings = stderr)
}

stanc_typed_ast <- function(stan_code, model_name = "model") {
  res <- run_stanc(stan_code, "--debug-decorated-ast", model_name)
  list(ast = paste(res$stdout, collapse = "\n"), warnings = res$warnings)
}

# inline #include files so line numbers refer to one self-contained program
stanc_expand_includes <- function(stan_code, include_paths, model_name = "model") {
  args <- c("--auto-format", "--canonicalize=includes")
  if (length(include_paths)) args <- c(args, paste0("--include-paths=", paste(include_paths, collapse = ",")))
  run_stanc(stan_code, args, model_name)$output  # the formatter writes to the -o file
}
