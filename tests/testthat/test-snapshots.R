# Generated R code for every committed model (or the error it raises), per mode.
# Review changes with testthat::snapshot_review().
models <- sort(list.files(test_path("models"), "\\.stan$", full.names = TRUE))

snapshot_text <- function(path, mode) {
  res <- tryCatch(suppressWarnings(stan2r(stan_program(path), mode, samps = if (mode != "prior_predictive") data.frame())),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
  res
}

for (f in models) {
  test_that(paste("translation of", basename(f)), {
    skip_if_no_stanc()
    expect_snapshot(cat(snapshot_text(f, "prior_predictive")))
  })
}
