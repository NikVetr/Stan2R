# Evaluate distribution templates from the mapping table as R functions.
dist_env <- function() {
  env <- new.env()
  for (h in Stan2R:::helper_source(names(Stan2R:::runtime_helpers))) eval(parse(text = h), envir = env)
  env
}

dist_fn <- function(name, which, args, env = dist_env()) {
  def <- Stan2R:::stan_distributions()[[name]]
  vals <- lapply(args, function(a) deparse(a))
  names(vals) <- def$args
  var <- switch(which, r = "n", d = "y", p = "y", q = "p")
  extra <- switch(which, r = list(n = "n"), d = list(y = "y"), p = list(y = "y", opts = ""), q = list(p = "p"))
  code <- Stan2R:::fill_template(def[[which]], c(vals, extra))
  f <- eval(parse(text = paste0("function(", var, ") ", code)))
  environment(f) <- env
  f
}

# typical arguments for each univariate distribution (Stan argument order)
dist_test_args <- list(
  normal = list(0.5, 2), std_normal = list(), student_t = list(4, 1, 2), cauchy = list(-1, 0.5),
  double_exponential = list(1, 2), logistic = list(1, 0.7), gumbel = list(1, 2),
  exp_mod_normal = list(0, 1, 2), skew_normal = list(0.5, 1.5, 3), lognormal = list(0.3, 0.6),
  chi_square = list(3), inv_chi_square = list(5), scaled_inv_chi_square = list(5, 0.8),
  exponential = list(1.5), gamma = list(2, 3), inv_gamma = list(3, 2), weibull = list(1.5, 2),
  frechet = list(3, 2), rayleigh = list(1.5), pareto = list(1, 3), pareto_type_2 = list(0, 2, 3),
  loglogistic = list(2, 4), beta = list(2, 5), beta_proportion = list(0.3, 10), uniform = list(-1, 2),
  bernoulli = list(0.3), bernoulli_logit = list(-0.5), binomial = list(10, 0.3), binomial_logit = list(10, 0.4),
  beta_binomial = list(12, 2, 3), hypergeometric = list(5, 7, 8), poisson = list(3.5), poisson_log = list(1.2),
  neg_binomial = list(3, 0.5), neg_binomial_2 = list(4, 2), neg_binomial_2_log = list(1.1, 3),
  discrete_range = list(2, 7)
)

# support used for numerical checks
dist_test_support <- function(name, def, args) {
  s <- def$support
  num <- function(v) if (is.na(suppressWarnings(as.numeric(v)))) args[[match(v, def$args)]] else as.numeric(v)
  c(num(s[1]), num(s[2]))
}
