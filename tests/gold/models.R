# Gold standard for whole models: the prior predictive from Stan2R's generated R
# code against Stan's own prior simulation. Each model has a hand-written
# tests/gold/models/<name>_prior.stan: the likelihood is removed, and observed
# outcomes are either sampled jointly as parameters (continuous outcomes) or
# drawn with _rng functions in generated quantities (discrete outcomes).
#
#   Rscript tests/gold/models.R [out_dir] [model ...]   (sequential compiles; 4 cores)

args <- commandArgs(TRUE)
out_dir <- if (length(args)) args[1] else "tmp/output/gold"
only <- args[-1]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(".", quiet = TRUE)
library(cmdstanr)
options(cmdstanr_no_ver_check = TRUE)

set.seed(1)
specs <- list(
  "beta-binomial" = list(
    data = list(n = 40, k = sample(5:60, 40, TRUE), n_groups = 6, group = sample(1:6, 40, TRUE)),
    rename = c(x = "x_sim"),
    qoi = function(g) list(theta_intercept = g("theta_intercept"), theta_group_sd = g("theta_group_sd"),
                           theta_group1 = g("theta_group")[, 1], conc = g("conc"), mu1 = g("mu")[, 1],
                           x1 = g("x")[, 1], x_mean = rowMeans(g("x")))),
  syn_regression_halft = list(
    data = list(N = 12, K = 3, X = matrix(round(rnorm(36), 2), 12, 3)),
    qoi = function(g) list(alpha = g("alpha"), beta2 = g("beta")[, 2], sigma = g("sigma"), mu1 = g("mu")[, 1],
                           y1 = g("y")[, 1], y_mean = rowMeans(g("y")), y_sd = apply(g("y"), 1, sd)),
    robust = c("sigma", "y1", "y_mean", "y_sd")),
  syn_hier_noncentered = list(
    data = list(J = 4, N = 10, g = c(1, 1, 2, 2, 3, 3, 4, 4, 1, 2)),
    qoi = function(g) list(mu = g("mu"), tau = g("tau"), theta1 = g("theta")[, 1], sigma = g("sigma"),
                           y1 = g("y")[, 1], y_mean = rowMeans(g("y")))),
  syn_truncation = list(
    data = list(N = 6, U = 3),
    qoi = function(g) list(mu = g("mu"), sigma = g("sigma"), y1 = g("y")[, 1], y_mean = rowMeans(g("y")))),
  syn_mixture_logmix = list(
    data = list(N = 8), rename = c(y = "y_sim"),
    qoi = function(g) list(lambda = g("lambda"), mu1 = g("mu")[, 1], mu2 = g("mu")[, 2], sigma = g("sigma"),
                           y1 = g("y")[, 1], y_mean = rowMeans(g("y")))),
  syn_mvn_lkj = list(
    data = list(N = 5, K = 3), rename = c(y = "y_sim"),
    qoi = function(g) list(mu1 = g("mu")[, 1], tau1 = g("tau")[, 1], omega12 = g("Omega")[, 1, 2],
                           omega23 = g("Omega")[, 2, 3], L22 = g("L")[, 2, 2], y11 = g("y")[, 1, 1],
                           y_sum = apply(g("y"), 1, sum)),
    robust = c("tau1", "y11", "y_sum")),
  syn_categorical = list(
    data = list(N = 20, K = 4, alpha = c(1, 2, 3, 4)), rename = c(z = "z_sim"),
    qoi = function(g) list(theta1 = g("theta")[, 1], theta4 = g("theta")[, 4], phi2 = g("phi")[, 2],
                           z1 = g("z")[, 1], z_mean = rowMeans(g("z")), z_new = g("z_new"))),
  syn_ordinal = list(
    data = list(N = 15, K = 4, x = round(seq(-2, 2, length.out = 15), 2)), rename = c(y = "y_sim"),
    qoi = function(g) list(beta = g("beta"), c1 = g("c")[, 1], c3 = g("c")[, 3], y1 = g("y")[, 1],
                           y_mean = rowMeans(g("y")))),
  syn_poisson_glm = list(
    data = list(N = 12, J = 3, g = rep(1:3, 4), log_exposure = round(log(seq(1, 4, length.out = 12)), 3)),
    rename = c(y = "y_sim"),
    qoi = function(g) list(alpha = g("alpha"), u1 = g("u")[, 1], u3 = g("u")[, 3], y1 = g("y")[, 1],
                           y_mean = rowMeans(g("y")))),
  student_t_mixture = list(
    data = NULL,
    qoi = function(g) list(x = g("x"), y = g("y")), robust = c("x", "y")),
  syn_functions = list(
    data = list(N = 6), rename = c(b = "b_sim"),
    qoi = function(g) list(a = g("a"), m = g("m"), y1 = g("y")[, 1], y_mean = rowMeans(g("y")),
                           b_mean = rowMeans(g("b")))),
  syn_reshape_matrix = list(
    data = list(R = 3, C = 2, w = c(0.5, -1)),
    qoi = function(g) list(Z12 = g("Z")[, 1, 2], row_effect2 = g("row_effect")[, 2], sigma = g("sigma"),
                           weighted1 = g("weighted")[, 1], grand = g("grand"), Y21 = g("Y")[, 2, 1])),
  syn_egtex_mixture = list(
    data = list(n = 10, total = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)), rename = c(count = "count_sim"),
    qoi = function(g) list(w1 = g("w")[, 1], conc = g("conc"), a = g("a"), count1 = g("count")[, 1],
                           count10 = g("count")[, 10], count_frac = rowMeans(g("count") / rep(c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50), each = nrow(g("count")))))),
  syn_affine = list(
    data = list(N = 5, nu = 6),
    qoi = function(g) list(r = g("r"), theta1 = g("theta")[, 1], theta_mean = rowMeans(g("theta")),
                           sigma = g("sigma"), s2_1 = g("s2")[, 1], s2_mean = rowMeans(g("s2")))),
  syn_spike_slab = list(
    data = list(n = 6, total = c(4, 8, 12, 16, 20, 40)), rename = c(count = "count_sim"),
    qoi = function(g) list(lambda = g("lambda"), theta1 = g("theta")[, 1], theta_mean = rowMeans(g("theta")),
                           count6 = g("count")[, 6], count_frac = rowMeans(g("count") / rep(c(4, 8, 12, 16, 20, 40), each = nrow(g("count")))))),
  syn_bounded_mixture = list(
    data = list(N = 4, trials = c(5, 10, 20, 40), s = 0.15, rate = 6), rename = c(successes = "successes_sim"),
    qoi = function(g) list(theta1 = g("theta")[, 1], theta_mean = rowMeans(g("theta")),
                           successes4 = g("successes")[, 4]))
)

# equal-probability Gamma precision nodes for the Student-t scale-mixture model
statmod_like_nodes <- function(K, df) {
  # precisions lambda ~ Gamma(df/2, df/2) discretised at equal-probability quantiles
  p <- (seq_len(K) - 0.5) / K
  list(lambda = qgamma(p, df / 2, df / 2), pi = rep(1 / K, K))
}
specs$student_t_mixture$data <- local({
  K <- 8; nodes <- statmod_like_nodes(K, 4)
  list(K = K, df = 4, loc = 1, scale = 2, lambda = nodes$lambda, pi = nodes$pi)
})

# draws-first arrays of each Stan variable
stan_getter <- function(fit, rename) {
  rv <- posterior::as_draws_rvars(fit$draws())
  function(name) {
    nm <- if (name %in% names(rename)) rename[[name]] else name
    x <- posterior::draws_of(rv[[nm]], with_chains = FALSE)
    if (length(dim(x)) == 2 && dim(x)[2] == 1) x <- x[, 1]
    x
  }
}
r_getter <- function(outs) {
  function(name) {
    vals <- lapply(outs, `[[`, name)
    first <- vals[[1]]
    if (length(first) == 1 && is.null(dim(first))) return(vapply(vals, as.numeric, 0))
    d <- if (is.null(dim(first))) length(first) else dim(first)
    arr <- array(unlist(vals), dim = c(d, length(vals)))
    aperm(arr, c(length(d) + 1, seq_along(d)))
  }
}

compare <- function(stan_x, stan_chains, r_x, name, robust) {
  probs <- c(0.05, 0.5, 0.95)
  sx <- posterior::as_draws_matrix(matrix(stan_x, ncol = 1, dimnames = list(NULL, name)))
  arr <- array(stan_x, dim = c(length(stan_x) / stan_chains, stan_chains))
  q_s <- quantile(stan_x, probs)
  q_r <- quantile(r_x, probs)
  se_q_s <- vapply(probs, function(p) posterior::mcse_quantile(arr, probs = p), 0)
  se_q_r <- vapply(probs, function(p) posterior::mcse_quantile(r_x, probs = p), 0)
  se <- sqrt(se_q_s^2 + se_q_r^2)
  # quantiles of discrete quantities can have zero Monte Carlo error; equal values then count as z = 0
  z <- ifelse(se > 0, (q_r - q_s) / se, ifelse(q_r == q_s, 0, NA))
  out <- data.frame(quantity = name, stat = c("q05", "q50", "q95"), stan = q_s, stan2r = q_r, z = z)
  if (!(name %in% robust)) {
    m <- data.frame(quantity = name, stat = c("mean", "sd"), stan = c(mean(stan_x), sd(stan_x)),
                    stan2r = c(mean(r_x), sd(r_x)))
    se_m <- c(posterior::mcse_mean(arr), posterior::mcse_sd(arr))
    se_r <- c(posterior::mcse_mean(r_x), posterior::mcse_sd(r_x))
    m$z <- (m$stan2r - m$stan) / sqrt(se_m^2 + se_r^2)
    out <- rbind(m, out)
  }
  out
}

results <- list()
diagnostics <- list()
for (name in names(specs)) {
  if (length(only) && !(name %in% only)) next
  sp <- specs[[name]]
  cat("==", name, "\n")
  model_file <- file.path("tests/testthat/models", paste0(name, ".stan"))
  gold_file <- file.path("tests/gold/models", paste0(name, "_prior.stan"))

  # Stan2R: prior predictive draws from the generated R code
  code <- suppressWarnings(parse_Stan(readLines(model_file), dat = sp$data, post_pred_sim = FALSE, sim = TRUE))
  writeLines(code, file.path(out_dir, paste0(name, "_prior_predictive.R")))
  expr <- parse(text = code)
  set.seed(42)
  n_r <- 8000
  outs <- vector("list", n_r)
  dat <- sp$data
  for (i in seq_len(n_r)) {
    env <- new.env()
    env$dat <- dat
    eval(expr, envir = env)
    outs[[i]] <- env$out
  }

  # Stan: prior simulation with the likelihood removed
  mod <- cmdstan_model(gold_file, quiet = TRUE, dir = out_dir)  # keep executables out of the repository
  stan_data <- sp$data[intersect(names(sp$data), names(mod$variables()$data))]
  fit <- mod$sample(data = stan_data, chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 2500,
                    adapt_delta = 0.95, seed = 7, refresh = 0, show_messages = FALSE, output_dir = out_dir)
  diag <- fit$diagnostic_summary(quiet = TRUE)
  summ <- fit$summary()
  diagnostics[[name]] <- data.frame(model = name, divergences = sum(diag$num_divergent),
                                    max_treedepth_hits = sum(diag$num_max_treedepth),
                                    max_rhat = max(summ$rhat, na.rm = TRUE),
                                    min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE))
  sg <- stan_getter(fit, sp$rename %||% character(0))
  rg <- r_getter(outs)
  q_stan <- sp$qoi(sg)
  q_r <- sp$qoi(rg)
  saveRDS(list(stan = q_stan, stan2r = q_r), file.path(out_dir, paste0(name, "_qoi.rds")))
  rows <- lapply(names(q_stan), function(q) compare(q_stan[[q]], 4, q_r[[q]], q, sp$robust %||% character(0)))
  res <- do.call(rbind, rows)
  res$model <- name
  results[[name]] <- res
  print(res[, c("quantity", "stat", "stan", "stan2r", "z")], row.names = FALSE, digits = 3)
}
res <- do.call(rbind, results)
diag <- do.call(rbind, diagnostics)
if (length(only) && file.exists(file.path(out_dir, "models.csv"))) {
  # a partial run replaces only its own models in the saved results
  old_res <- read.csv(file.path(out_dir, "models.csv"))
  old_diag <- read.csv(file.path(out_dir, "models_diagnostics.csv"))
  res <- rbind(old_res[!(old_res$model %in% names(results)), names(res)], res)
  diag <- rbind(old_diag[!(old_diag$model %in% names(results)), names(diag)], diag)
}
write.csv(res, file.path(out_dir, "models.csv"), row.names = FALSE)
write.csv(diag, file.path(out_dir, "models_diagnostics.csv"), row.names = FALSE)
cat("\nStan diagnostics:\n")
print(diag, row.names = FALSE)
by_model <- do.call(rbind, lapply(split(res, res$model), function(r)
  data.frame(model = r$model[1], comparisons = nrow(r), max_abs_z = max(abs(r$z), na.rm = TRUE),
             n_abs_z_gt_3 = sum(abs(r$z) > 3, na.rm = TRUE), undefined_z = sum(is.na(r$z)))))
write.csv(by_model, file.path(out_dir, "models_summary.csv"), row.names = FALSE)
print(by_model, row.names = FALSE)
cat(sprintf("\n%d comparisons over %d models; max |z| = %.2f; |z| > 3: %d; |z| > 4: %d; undefined: %d\n",
            nrow(res), length(unique(res$model)), max(abs(res$z), na.rm = TRUE), sum(abs(res$z) > 3, na.rm = TRUE),
            sum(abs(res$z) > 4, na.rm = TRUE), sum(is.na(res$z))))
