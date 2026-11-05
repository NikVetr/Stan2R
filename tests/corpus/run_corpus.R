# Convert every .stan file in the given directories with Stan2R and run the
# generated R code on mock data from stanc (--debug-generate-data/-inits).
#
#   Rscript tests/corpus/run_corpus.R --out tmp/output/corpus <dir-or-file> ...
#
# Writes results.csv (one row per model and mode) and summary.md.

args <- commandArgs(TRUE)
out_i <- which(args == "--out")
out_dir <- if (length(out_i)) args[out_i + 1] else "tmp/output/corpus"
inputs <- if (length(out_i)) args[-c(out_i, out_i + 1)] else args
if (!length(inputs)) stop("give at least one model directory or .stan file")
pkgload::load_all(".", quiet = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- unlist(lapply(inputs, function(p) if (dir.exists(p)) list.files(p, "\\.stan$", full.names = TRUE) else p))
modes <- c("prior_predictive", "posterior_predictive", "posterior_density")

stanc_json <- function(path, flag, data_file = NULL) {
  out <- tempfile(fileext = ".json")
  extra <- if (!is.null(data_file)) paste0("--debug-data-file=", data_file) else character(0)
  status <- system2(Stan2R:::stanc_path(), c(flag, extra, "-o", out, shQuote(path)), stdout = FALSE, stderr = FALSE)
  if (status != 0 || !file.exists(out)) return(NULL)
  attr_file <- out
  structure(jsonlite::fromJSON(out, simplifyVector = TRUE), file = attr_file)
}

# one CmdStan-style draws row from a list of values
draw_row <- function(values) {
  cols <- list()
  for (nm in names(values)) {
    v <- values[[nm]]
    if (length(dim(v)) >= 2 || length(v) > 1 || is.list(v)) {
      v <- if (is.list(v)) do.call(rbind, v) else v
      if (is.null(dim(v))) dim(v) <- length(v)
      idx <- arrayInd(seq_along(v), dim(v))
      names_v <- paste0(nm, "[", apply(idx, 1, paste, collapse = ","), "]")
      cols[names_v] <- as.list(as.vector(v))
    } else {
      cols[[nm]] <- as.numeric(v)
    }
  }
  as.data.frame(cols, check.names = FALSE)
}

classify <- function(e) {
  if (inherits(e, "stan2r_unsupported")) return(c(status = "unsupported", construct = e$construct))
  msg <- conditionMessage(e)
  if (grepl("stanc rejected", msg)) return(c(status = "stanc_error", construct = "rejected by stanc 2.38"))
  c(status = "error", construct = substr(gsub("\\s+", " ", msg), 1, 120))
}

rows <- list()
for (f in files) {
  id <- sub("\\.stan$", "", basename(f))
  program <- tryCatch(stan_program(f), error = function(e) e)
  dat <- stanc_json(f, "--debug-generate-data")
  inits <- if (!is.null(dat)) stanc_json(f, "--debug-generate-inits", attr(dat, "file")) else NULL
  prior_out <- NULL
  for (mode in modes) {
    row <- list(model = id, mode = mode, status = "ok", construct = "", runs = NA, lines = NA_integer_, warnings = "")
    if (inherits(program, "error")) {
      cl <- classify(program)
      row[names(cl)] <- as.list(cl)
      rows[[length(rows) + 1]] <- row
      next
    }
    samps <- if (!is.null(inits)) tryCatch(draw_row(inits), error = function(e) NULL) else NULL
    if (is.null(samps) && !is.null(prior_out)) {
      # no mock inits from stanc: use the parameters of the prior predictive run
      params <- program$symbols$name[program$symbols$role == "param"]
      samps <- tryCatch(draw_row(prior_out[params]), error = function(e) NULL)
    }
    code <- tryCatch(withCallingHandlers(
      stan2r(program, mode, dat = dat, samps = if (mode == "prior_predictive") NULL else samps),
      warning = function(w) invokeRestart("muffleWarning")), error = function(e) e)
    if (inherits(code, "error")) {
      cl <- classify(code)
      row[names(cl)] <- as.list(cl)
    } else {
      row$lines <- length(strsplit(code, "\n")[[1]])
      row$warnings <- paste(attr(code, "warnings"), collapse = " | ")
      writeLines(code, file.path(out_dir, paste0(id, "__", mode, ".R")))
      if (!is.null(dat)) {
        env <- new.env()
        env$dat <- dat
        set.seed(1)
        run <- tryCatch({ suppressWarnings(eval(parse(text = code), envir = env)); TRUE },
                        error = function(e) conditionMessage(e))
        row$runs <- isTRUE(run)
        if (isTRUE(run) && mode == "prior_predictive") prior_out <- env$out
        if (!isTRUE(run)) row$construct <- paste("runtime:", substr(gsub("\\s+", " ", run), 1, 100))
      }
    }
    rows[[length(rows) + 1]] <- row
  }
}
res <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
write.csv(res, file.path(out_dir, "results.csv"), row.names = FALSE)

md <- c("# Stan2R corpus run", "", sprintf("%d models, stanc %s, %s", length(files), Stan2R:::stanc_version_string(), Sys.Date()), "")
for (m in modes) {
  r <- res[res$mode == m, ]
  md <- c(md, sprintf("- %s: %d translated (%d of those ran on mock data), %d unsupported, %d stanc errors, %d other errors",
                      m, sum(r$status == "ok"), sum(r$runs %in% TRUE), sum(r$status == "unsupported"),
                      sum(r$status == "stanc_error"), sum(r$status == "error")))
}
md <- c(md, "", "## Failures by construct (prior_predictive)", "", "| construct | models |", "| --- | --- |")
r <- res[res$mode == "prior_predictive" & (res$status != "ok" | res$runs %in% FALSE), ]
for (k in unique(r$construct)) md <- c(md, sprintf("| %s | %s |", k, paste(r$model[r$construct == k], collapse = ", ")))
writeLines(md, file.path(out_dir, "summary.md"))
cat(md, sep = "\n")
