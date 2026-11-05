test_that("types, bounds and transforms come from stanc", {
  skip_if_no_stanc()
  p <- stan_program(paste(
    "data { int<lower=0> N; matrix[N, 2] X; array[N] vector[3] v; }",
    "parameters { real<lower=0, upper=1> p; simplex[3] s; cholesky_factor_corr[2] L; }",
    "model { p ~ beta(2, 2); }", sep = "\n"))
  sy <- p$symbols
  expect_equal(sy$name, c("N", "X", "v", "p", "s", "L"))
  expect_equal(sy$base, c("int", "matrix", "vector", "real", "vector", "matrix"))
  expect_equal(sy$rank, c(0L, 2L, 2L, 0L, 1L, 2L))
  expect_equal(sy$transform, c("bounds", "none", "none", "bounds", "simplex", "cholesky_factor_corr"))
})

test_that("statements get original line numbers and comments", {
  skip_if_no_stanc()
  code <- c("data {",
            "  int N; vector[N] y; // two on a line",
            "}",
            "parameters {",
            "  /* the mean",
            "     of y */",
            "  real mu;",
            "}",
            "model {",
            "  for (n in 1:N)",
            "    y[n] ~ normal(mu,",
            "                  1);  // likelihood",
            "  if (N > 1) mu ~ normal(0, 1); else mu ~ normal(0, 2);",
            "}")
  p <- stan_program(code)
  d <- p$blocks$data
  expect_equal(vapply(d, `[[`, 0L, "line"), c(2L, 2L))
  expect_equal(d[[2]]$comment_after, "two on a line")
  expect_equal(p$blocks$parameters[[1]]$comments_before[[1]]$text, c("the mean", "of y"))
  loop <- p$blocks$model[[1]]
  expect_equal(c(loop$line, loop$end_line), c(10L, 12L))
  expect_equal(loop$body$line, 11L)
  expect_equal(loop$body$comment_after, "likelihood")
  branch <- p$blocks$model[[2]]
  expect_equal(branch$node, "If")
  expect_equal(branch$otherwise$line, 13L)
})

test_that("stanc errors are reported with stanc's message", {
  skip_if_no_stanc()
  expect_error(stan_program("parameters { real x }"), "stanc rejected the model")
})

test_that("#include files are expanded before parsing", {
  skip_if_no_stanc()
  dir <- withr::local_tempdir()
  writeLines("real twice(real x) { return 2 * x; }", file.path(dir, "f.stanfunctions"))
  p <- stan_program(c("functions {", "#include f.stanfunctions", "}", "parameters { real x; }",
                      "model { x ~ normal(twice(1), 1); }"), include_paths = dir)
  expect_true("twice" %in% names(p$functions))
})
