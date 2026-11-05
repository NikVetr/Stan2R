bb_data <- function(n = 50, n_groups = 5, seed = 1) {
  set.seed(seed)
  list(n = n, k = sample(1:30, n, replace = TRUE), x = rep(0L, n), n_groups = n_groups,
       group = sample(1:n_groups, n, replace = TRUE))
}

test_that("the README beta-binomial model runs in every mode", {
  skip_if_no_stanc()
  stan_code <- readLines(model_path("beta-binomial"))
  dat <- bb_data()
  prior <- run_code(parse_Stan(stan_code, dat = dat, post_pred_sim = FALSE, sim = TRUE), dat)
  expect_named(prior, c("n", "k", "x", "n_groups", "group", "k_scaled", "theta_intercept", "theta_slope",
                        "theta_group", "theta_group_sd", "conc_offset", "theta", "mu", "conc", "shape1", "shape2"))
  expect_true(all(prior$x >= 0 & prior$x <= dat$k))
  expect_gt(prior$theta_group_sd, 0)

  samps <- data.frame(theta_intercept = 0.2, theta_slope = -0.1, theta_group_sd = 0.5, conc_offset = 2)
  for (j in 1:5) samps[[paste0("theta_group.", j)]] <- j / 10
  post <- run_code(parse_Stan(stan_code, dat = dat, samps = samps, sample_index = 1), dat)
  expect_equal(post$theta_group, (1:5) / 10)
  dens <- run_code(parse_Stan(stan_code, dat = dat, samps = samps, sample_index = 1, sim = FALSE), dat)
  shape1 <- dens$shape1
  shape2 <- dens$shape2
  expect_equal(dens$x, choose(dat$k, dat$x) * beta(dat$x + shape1, dat$k - dat$x + shape2) / beta(shape1, shape2))
})

test_that("CmdStan column names with brackets are read too", {
  skip_if_no_stanc()
  code <- "parameters { matrix[2, 2] m; vector[3] v; } model { to_vector(m) ~ std_normal(); v ~ std_normal(); }"
  samps <- data.frame(`m[1,1]` = 1, `m[2,1]` = 2, `m[1,2]` = 3, `m[2,2]` = 4, `v[1]` = 5, `v[2]` = 6, `v[3]` = 7,
                      check.names = FALSE)
  out <- run_code(tr(code, "posterior_predictive", samps = samps), list())
  expect_equal(out$m, matrix(1:4, 2))
  expect_equal(out$v, 5:7)
})

test_that("bounded priors are truncated, and student_t keeps location and scale", {
  skip_if_no_stanc()
  code <- "parameters { real<lower=0> s; real<lower=1, upper=3> u; real t; }
           model { s ~ student_t(3, 0, 2.5); u ~ normal(0, 5); t ~ student_t(5, 10, 0.1); }"
  lines <- code_lines(tr(code))
  expect_true("s <- abs(0 + 2.5 * rt(1, 3))" %in% lines)
  expect_true("u <- rtruncated(1, function(q, ...) pnorm(q, 0, 5, ...), function(p, ...) qnorm(p, 0, 5, ...), 1, 3)" %in% lines)
  draws <- replicate(2000, unlist(run_code(tr(code), list(), seed = NULL)))
  expect_true(all(draws["s", ] > 0))
  expect_true(all(draws["u", ] > 1 & draws["u", ] < 3))
  expect_equal(median(draws["t", ]), 10, tolerance = 0.01)
  expect_equal(median(draws["s", ]), 2.5 * qt(0.75, 3), tolerance = 0.1)
})

test_that("a truncation whose normaliser depends on parameters is flagged", {
  skip_if_no_stanc()
  code <- "parameters { real m; real<lower=0> theta; } model { m ~ normal(0, 1); theta ~ normal(m, 1); }"
  expect_warning(stan2r(stan_program(code), "prior_predictive"), "normalising constant depends on other parameters")
  code_ok <- "parameters { real<lower=0> tau; real<lower=0> s; } model { tau ~ normal(0, 1); s ~ normal(0, tau); }"
  expect_no_warning(stan2r(stan_program(code_ok), "prior_predictive"))
})

test_that("prior draws are scheduled before the quantities that use them", {
  skip_if_no_stanc()
  code <- "parameters { real log_tau; vector[3] theta; }
           transformed parameters { real tau = exp(log_tau); }
           model { theta ~ normal(0, tau); log_tau ~ normal(0, 1); }"
  lines <- code_lines(tr(code))
  pos <- function(pattern) grep(pattern, lines)[1]
  expect_lt(pos("^log_tau <- rnorm"), pos("^tau <- exp"))
  expect_lt(pos("^tau <- exp"), pos("^theta <- rnorm"))
})

test_that("mixtures draw a component, and density mode mixes densities", {
  skip_if_no_stanc()
  code <- "data { int N; vector[N] y; } parameters { real<lower=0, upper=1> w; }
           model { for (n in 1:N) target += log_sum_exp(log(w) + normal_lpdf(y[n] | -5, 1), log1m(w) + normal_lpdf(y[n] | 5, 1)); }"
  dat <- list(N = 2000, y = rep(0, 2000))
  out <- run_code(tr(code, dat = dat), dat)
  expect_equal(mean(out$y > 0), 1 - out$w, tolerance = 0.05)
  samps <- data.frame(w = 0.25)
  dens <- run_code(tr(code, "posterior_density", dat = dat, samps = samps), dat)
  expect_equal(dens$y[1], 0.25 * dnorm(0, -5) + 0.75 * dnorm(0, 5))
})

test_that("unsupported constructs stop with the construct and its line", {
  skip_if_no_stanc()
  cases <- list(
    list(code = "data { vector[2] y; } parameters { real mu; } model { y ~ normal(mu, 1); }",
         construct = "parameter without a prior", line = 1L),
    list(code = "parameters { real mu; }\ntransformed parameters { real e = exp(mu); }\nmodel {\n  e ~ normal(0, 1);\n}",
         construct = "~ on transformed parameter", line = 4L),
    list(code = "parameters { real mu; }\nmodel {\n  mu ~ normal(0, 1);\n  target += -mu^4;\n}",
         construct = "free-form target increment", line = 4L),
    list(code = "parameters { real<lower=0> s; }\nmodel {\n  log(s) ~ normal(0, 1);\n}",
         construct = "expression on left of ~", line = 3L))
  for (cs in cases) {
    err <- tryCatch(tr(cs$code), error = function(e) e)
    expect_s3_class(err, "stan2r_unsupported")
    expect_equal(err$construct, cs$construct)
    expect_equal(err$line, cs$line)
  }
})

test_that("posterior modes skip priors and keep transformed parameters", {
  skip_if_no_stanc()
  code <- "parameters { real log_s; } transformed parameters { real s = exp(log_s); }
           model { target += -log_s^2; }"
  out <- run_code(tr(code, "posterior_predictive", samps = data.frame(log_s = log(2))), list())
  expect_equal(out$s, 2)
})

test_that("user-defined functions, _rng and derived draws", {
  skip_if_no_stanc()
  stan_code <- readLines(model_path("syn_functions"))
  dat <- list(N = 4000, y = rep(0, 4000), b = rep(0L, 4000))
  code <- parse_Stan(stan_code, dat = dat, post_pred_sim = FALSE)
  expect_match(code, "two_normals_draw__ <- function(y, w, m)", fixed = TRUE)
  expect_match(code, "m <- rtruncated_reject(1, function(n) shifted_exp_rng(0.5, 2), 0, Inf)", fixed = TRUE)
  out <- run_code(code, dat)
  expect_gt(out$m, 0.5)
  expect_equal(mean(out$y < 0), 0.3 * pnorm(out$m) + 0.7 * pnorm(-out$m), tolerance = 0.05)
})

test_that("affine left-hand sides are inverted, with Jacobians and bounds", {
  skip_if_no_stanc()
  stan_code <- readLines(model_path("syn_affine"))
  dat <- list(N = 3, s2 = c(1, 2, 3), nu = 4)
  expect_warning(code <- stan2r(stan_program(stan_code), "prior_predictive", dat = dat), "parameter-dependent factor")
  draws <- replicate(3000, unlist(run_code(code, dat, seed = NULL)[c("r", "theta", "sigma")]))
  expect_true(all(abs(draws["r", ]) < 1))
  expect_true(all(draws[2:4, ] > 0 & draws[2:4, ] < 1))
  expect_equal(unname(quantile(draws["r", ], c(0.1, 0.5, 0.9))), 2 * qbeta(c(0.1, 0.5, 0.9), 2, 2) - 1, tolerance = 0.05)
  dens_code <- suppressWarnings(stan2r(stan_program(stan_code), "posterior_density", dat = dat,
                                       samps = data.frame(r = 0, `theta[1]` = 0.5, `theta[2]` = 0.5, `theta[3]` = 0.5,
                                                          sigma = 1.5, check.names = FALSE)))
  out <- run_code(dens_code, dat)
  x <- 4 * dat$s2 / 1.5^2
  expect_equal(out$s2, dchisq(x, 4) * 4 / 1.5^2)
})

test_that("likelihood terms reached through local variables are not dropped", {
  skip_if_no_stanc()
  code <- "data { int N; vector[N] y; } parameters { real mu; }
           model { mu ~ normal(0, 1); { vector[N] ll; for (n in 1:N) ll[n] = -0.5 * square(y[n] - mu); target += sum(ll); } }"
  err <- tryCatch(tr(code, "posterior_density", samps = data.frame(mu = 0)), error = function(e) e)
  expect_s3_class(err, "stan2r_unsupported")
  expect_equal(err$construct, "free-form target increment")
  code_tp <- "data { int N; vector[N] y; } parameters { real mu; } transformed parameters { vector[N] r = y - mu; }
              model { mu ~ normal(0, 1); r ~ normal(0, 1); }"
  err <- tryCatch(tr(code_tp, "posterior_predictive", samps = data.frame(mu = 0)), error = function(e) e)
  expect_s3_class(err, "stan2r_unsupported")
  expect_match(conditionMessage(err), "computed from observed data")
})
