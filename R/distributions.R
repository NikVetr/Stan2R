# Stan distributions -> R code templates, following the parameterisations in the
# Stan Functions Reference. Placeholders: {y} outcome, {n} number of draws,
# {p} probability (quantile functions), {opts} extra arguments for native R CDFs
# (lower.tail / log.p), and each Stan argument by name. `{(name)}` inserts the
# argument in parentheses when it is not atomic, for use inside arithmetic.
#
#   r: random draws         d: log density/mass (elementwise, log = TRUE)
#   p: CDF (elementwise)    q: quantile (inverse CDF), used for truncation
#   p_native: TRUE when `p` is an R function that accepts lower.tail / log.p
#   support: c(lower, upper); numbers, or a Stan argument name
#   helpers: runtime helper functions the templates call (see runtime.R)
#
# Multivariate entries have `multivariate = TRUE`; their r templates receive the
# outcome as {like} (an object of the outcome's shape to fill).

dist <- function(args, r = NULL, d = NULL, p = NULL, q = NULL, p_native = FALSE,
                 support = c(-Inf, Inf), discrete = FALSE, multivariate = FALSE,
                 helpers = character(0), note = NULL) {
  list(args = args, r = r, d = d, p = p, q = q, p_native = p_native, support = support,
       discrete = discrete, multivariate = multivariate, helpers = helpers, note = note)
}

stan_distributions <- function() {
  list(
    # unbounded continuous
    normal = dist(c("mu", "sigma"), r = "rnorm({n}, {mu}, {sigma})",
                  d = "dnorm({y}, {mu}, {sigma}, log = TRUE)", p = "pnorm({y}, {mu}, {sigma}{opts})",
                  q = "qnorm({p}, {mu}, {sigma}{opts})", p_native = TRUE),
    std_normal = dist(character(0), r = "rnorm({n})", d = "dnorm({y}, log = TRUE)",
                      p = "pnorm({y}{opts})", q = "qnorm({p}{opts})", p_native = TRUE),
    student_t = dist(c("nu", "mu", "sigma"), r = "{(mu)} + {(sigma)} * rt({n}, {nu})",
                     d = "dt(({y} - {(mu)}) / {(sigma)}, {nu}, log = TRUE) - log({sigma})",
                     p = "pt(({y} - {(mu)}) / {(sigma)}, {nu}{opts})",
                     q = "{(mu)} + {(sigma)} * qt({p}, {nu}{opts})", p_native = TRUE),
    cauchy = dist(c("mu", "sigma"), r = "rcauchy({n}, {mu}, {sigma})",
                  d = "dcauchy({y}, {mu}, {sigma}, log = TRUE)", p = "pcauchy({y}, {mu}, {sigma}{opts})",
                  q = "qcauchy({p}, {mu}, {sigma}{opts})", p_native = TRUE),
    double_exponential = dist(c("mu", "sigma"), r = "{(mu)} + {(sigma)} * (rexp({n}) - rexp({n}))",
                              d = "-log(2 * {(sigma)}) - abs({y} - {(mu)}) / {(sigma)}",
                              p = "plaplace({y}, {mu}, {sigma})", q = "qlaplace({p}, {mu}, {sigma})",
                              helpers = c("plaplace", "qlaplace")),
    logistic = dist(c("mu", "sigma"), r = "rlogis({n}, {mu}, {sigma})",
                    d = "dlogis({y}, {mu}, {sigma}, log = TRUE)", p = "plogis({y}, {mu}, {sigma}{opts})",
                    q = "qlogis({p}, {mu}, {sigma}{opts})", p_native = TRUE),
    gumbel = dist(c("mu", "beta"), r = "{(mu)} - {(beta)} * log(-log(runif({n})))",
                  d = "-log({beta}) - ({y} - {(mu)}) / {(beta)} - exp(-({y} - {(mu)}) / {(beta)})",
                  p = "exp(-exp(-({y} - {(mu)}) / {(beta)}))", q = "{(mu)} - {(beta)} * log(-log({p}))"),
    skew_double_exponential = dist(c("mu", "sigma", "tau"),
                                   r = "qskew_double_exponential(runif({n}), {mu}, {sigma}, {tau})",
                                   d = "log(2) + log({tau}) + log1p(-{(tau)}) - log({sigma}) - 2 / {(sigma)} * ifelse({y} < {(mu)}, (1 - {(tau)}) * ({(mu)} - {y}), {(tau)} * ({y} - {(mu)}))",
                                   p = "pskew_double_exponential({y}, {mu}, {sigma}, {tau})",
                                   q = "qskew_double_exponential({p}, {mu}, {sigma}, {tau})",
                                   helpers = c("pskew_double_exponential", "qskew_double_exponential")),
    exp_mod_normal = dist(c("mu", "sigma", "lambda"), r = "rnorm({n}, {mu}, {sigma}) + rexp({n}, {lambda})",
                          d = "dexp_mod_normal({y}, {mu}, {sigma}, {lambda})",
                          p = "pexp_mod_normal({y}, {mu}, {sigma}, {lambda})",
                          helpers = c("dexp_mod_normal", "pexp_mod_normal")),
    skew_normal = dist(c("xi", "omega", "alpha"), r = "rskew_normal({n}, {xi}, {omega}, {alpha})",
                       d = "log(2) - log({omega}) + dnorm(({y} - {(xi)}) / {(omega)}, log = TRUE) + pnorm({(alpha)} * ({y} - {(xi)}) / {(omega)}, log.p = TRUE)",
                       helpers = "rskew_normal"),
    # positive continuous
    lognormal = dist(c("mu", "sigma"), r = "rlnorm({n}, {mu}, {sigma})",
                     d = "dlnorm({y}, {mu}, {sigma}, log = TRUE)", p = "plnorm({y}, {mu}, {sigma}{opts})",
                     q = "qlnorm({p}, {mu}, {sigma}{opts})", p_native = TRUE, support = c(0, Inf)),
    chi_square = dist("nu", r = "rchisq({n}, {nu})", d = "dchisq({y}, {nu}, log = TRUE)",
                      p = "pchisq({y}, {nu}{opts})", q = "qchisq({p}, {nu}{opts})", p_native = TRUE,
                      support = c(0, Inf)),
    inv_chi_square = dist("nu", r = "1 / rchisq({n}, {nu})",
                          d = "dchisq(1 / {y}, {nu}, log = TRUE) - 2 * log({y})",
                          p = "pchisq(1 / {y}, {nu}, lower.tail = FALSE)", q = "1 / qchisq(1 - {p}, {nu})",
                          support = c(0, Inf)),
    scaled_inv_chi_square = dist(c("nu", "s"), r = "{(nu)} * {(s)}^2 / rchisq({n}, {nu})",
                                 d = "dchisq({(nu)} * {(s)}^2 / {y}, {nu}, log = TRUE) + log({(nu)} * {(s)}^2) - 2 * log({y})",
                                 p = "pchisq({(nu)} * {(s)}^2 / {y}, {nu}, lower.tail = FALSE)",
                                 q = "{(nu)} * {(s)}^2 / qchisq(1 - {p}, {nu})", support = c(0, Inf)),
    exponential = dist("beta", r = "rexp({n}, {beta})", d = "dexp({y}, {beta}, log = TRUE)",
                       p = "pexp({y}, {beta}{opts})", q = "qexp({p}, {beta}{opts})", p_native = TRUE,
                       support = c(0, Inf)),
    gamma = dist(c("alpha", "beta"), r = "rgamma({n}, shape = {alpha}, rate = {beta})",
                 d = "dgamma({y}, shape = {alpha}, rate = {beta}, log = TRUE)",
                 p = "pgamma({y}, shape = {alpha}, rate = {beta}{opts})",
                 q = "qgamma({p}, shape = {alpha}, rate = {beta}{opts})", p_native = TRUE, support = c(0, Inf)),
    inv_gamma = dist(c("alpha", "beta"), r = "1 / rgamma({n}, shape = {alpha}, rate = {beta})",
                     d = "dgamma(1 / {y}, shape = {alpha}, rate = {beta}, log = TRUE) - 2 * log({y})",
                     p = "pgamma(1 / {y}, shape = {alpha}, rate = {beta}, lower.tail = FALSE)",
                     q = "1 / qgamma(1 - {p}, shape = {alpha}, rate = {beta})", support = c(0, Inf)),
    weibull = dist(c("alpha", "sigma"), r = "rweibull({n}, shape = {alpha}, scale = {sigma})",
                   d = "dweibull({y}, shape = {alpha}, scale = {sigma}, log = TRUE)",
                   p = "pweibull({y}, shape = {alpha}, scale = {sigma}{opts})",
                   q = "qweibull({p}, shape = {alpha}, scale = {sigma}{opts})", p_native = TRUE, support = c(0, Inf)),
    frechet = dist(c("alpha", "sigma"), r = "1 / rweibull({n}, shape = {alpha}, scale = 1 / {(sigma)})",
                   d = "log({alpha}) - log({sigma}) - (1 + {(alpha)}) * log({y} / {(sigma)}) - ({y} / {(sigma)})^(-{(alpha)})",
                   p = "exp(-({y} / {(sigma)})^(-{(alpha)}))", q = "{(sigma)} * (-log({p}))^(-1 / {(alpha)})",
                   support = c(0, Inf)),
    rayleigh = dist("sigma", r = "{(sigma)} * sqrt(-2 * log(runif({n})))",
                    d = "log({y}) - 2 * log({sigma}) - {y}^2 / (2 * {(sigma)}^2)",
                    p = "1 - exp(-{y}^2 / (2 * {(sigma)}^2))", q = "{(sigma)} * sqrt(-2 * log(1 - {p}))",
                    support = c(0, Inf)),
    pareto = dist(c("y_min", "alpha"), r = "{(y_min)} * runif({n})^(-1 / {(alpha)})",
                  d = "ifelse({y} < {y_min}, -Inf, log({alpha}) + {(alpha)} * log({y_min}) - ({(alpha)} + 1) * log({y}))",
                  p = "pmax(0, 1 - ({(y_min)} / {y})^{(alpha)})", q = "{(y_min)} * (1 - {p})^(-1 / {(alpha)})",
                  support = c("y_min", Inf)),
    pareto_type_2 = dist(c("mu", "lambda", "alpha"),
                         r = "{(mu)} + {(lambda)} * (runif({n})^(-1 / {(alpha)}) - 1)",
                         d = "log({alpha}) - log({lambda}) - ({(alpha)} + 1) * log1p(({y} - {(mu)}) / {(lambda)})",
                         p = "1 - (1 + ({y} - {(mu)}) / {(lambda)})^(-{(alpha)})",
                         q = "{(mu)} + {(lambda)} * ((1 - {p})^(-1 / {(alpha)}) - 1)", support = c("mu", Inf)),
    loglogistic = dist(c("alpha", "beta"), r = "{(alpha)} * (1 / runif({n}) - 1)^(-1 / {(beta)})",
                       d = "log({beta}) - log({alpha}) + ({(beta)} - 1) * log({y} / {(alpha)}) - 2 * log1p(({y} / {(alpha)})^{(beta)})",
                       p = "1 / (1 + ({y} / {(alpha)})^(-{(beta)}))",
                       q = "{(alpha)} * ({p} / (1 - {p}))^(1 / {(beta)})", support = c(0, Inf)),
    # bounded continuous
    beta = dist(c("alpha", "beta"), r = "rbeta({n}, {alpha}, {beta})",
                d = "dbeta({y}, {alpha}, {beta}, log = TRUE)", p = "pbeta({y}, {alpha}, {beta}{opts})",
                q = "qbeta({p}, {alpha}, {beta}{opts})", p_native = TRUE, support = c(0, 1)),
    beta_proportion = dist(c("mu", "kappa"), r = "rbeta({n}, {(mu)} * {(kappa)}, (1 - {(mu)}) * {(kappa)})",
                           d = "dbeta({y}, {(mu)} * {(kappa)}, (1 - {(mu)}) * {(kappa)}, log = TRUE)",
                           p = "pbeta({y}, {(mu)} * {(kappa)}, (1 - {(mu)}) * {(kappa)}{opts})",
                           q = "qbeta({p}, {(mu)} * {(kappa)}, (1 - {(mu)}) * {(kappa)}{opts})", p_native = TRUE,
                           support = c(0, 1)),
    uniform = dist(c("alpha", "beta"), r = "runif({n}, {alpha}, {beta})",
                   d = "dunif({y}, {alpha}, {beta}, log = TRUE)", p = "punif({y}, {alpha}, {beta}{opts})",
                   q = "qunif({p}, {alpha}, {beta}{opts})", p_native = TRUE, support = c("alpha", "beta")),
    # discrete
    bernoulli = dist("theta", r = "rbinom({n}, 1, {theta})", d = "dbinom({y}, 1, {theta}, log = TRUE)",
                     p = "pbinom({y}, 1, {theta}{opts})", q = "qbinom({p}, 1, {theta}{opts})", p_native = TRUE,
                     discrete = TRUE, support = c(0, 1)),
    bernoulli_logit = dist("alpha", r = "rbinom({n}, 1, plogis({alpha}))",
                           d = "dbinom({y}, 1, plogis({alpha}), log = TRUE)",
                           p = "pbinom({y}, 1, plogis({alpha}){opts})", q = "qbinom({p}, 1, plogis({alpha}){opts})",
                           p_native = TRUE, discrete = TRUE, support = c(0, 1)),
    binomial = dist(c("N", "theta"), r = "rbinom({n}, {N}, {theta})", d = "dbinom({y}, {N}, {theta}, log = TRUE)",
                    p = "pbinom({y}, {N}, {theta}{opts})", q = "qbinom({p}, {N}, {theta}{opts})", p_native = TRUE,
                    discrete = TRUE, support = c(0, "N")),
    binomial_logit = dist(c("N", "alpha"), r = "rbinom({n}, {N}, plogis({alpha}))",
                          d = "dbinom({y}, {N}, plogis({alpha}), log = TRUE)",
                          p = "pbinom({y}, {N}, plogis({alpha}){opts})", q = "qbinom({p}, {N}, plogis({alpha}){opts})",
                          p_native = TRUE, discrete = TRUE, support = c(0, "N")),
    beta_binomial = dist(c("N", "alpha", "beta"), r = "rbinom({n}, {N}, rbeta({n}, {alpha}, {beta}))",
                         d = "lchoose({N}, {y}) + lbeta({y} + {(alpha)}, {(N)} - {y} + {(beta)}) - lbeta({alpha}, {beta})",
                         p = "pbeta_binomial({y}, {N}, {alpha}, {beta})", q = "qbeta_binomial({p}, {N}, {alpha}, {beta})",
                         discrete = TRUE, support = c(0, "N"), helpers = c("pbeta_binomial", "qbeta_binomial")),
    hypergeometric = dist(c("N", "a", "b"), r = "rhyper({n}, {a}, {b}, {N})",
                          d = "dhyper({y}, {a}, {b}, {N}, log = TRUE)", p = "phyper({y}, {a}, {b}, {N}{opts})",
                          q = "qhyper({p}, {a}, {b}, {N}{opts})", p_native = TRUE, discrete = TRUE, support = c(0, "N")),
    poisson = dist("lambda", r = "rpois({n}, {lambda})", d = "dpois({y}, {lambda}, log = TRUE)",
                   p = "ppois({y}, {lambda}{opts})", q = "qpois({p}, {lambda}{opts})", p_native = TRUE,
                   discrete = TRUE, support = c(0, Inf)),
    poisson_log = dist("alpha", r = "rpois({n}, exp({alpha}))", d = "dpois({y}, exp({alpha}), log = TRUE)",
                       p = "ppois({y}, exp({alpha}){opts})", q = "qpois({p}, exp({alpha}){opts})", p_native = TRUE,
                       discrete = TRUE, support = c(0, Inf)),
    neg_binomial = dist(c("alpha", "beta"), r = "rnbinom({n}, size = {alpha}, prob = {(beta)} / (1 + {(beta)}))",
                        d = "dnbinom({y}, size = {alpha}, prob = {(beta)} / (1 + {(beta)}), log = TRUE)",
                        p = "pnbinom({y}, size = {alpha}, prob = {(beta)} / (1 + {(beta)}){opts})",
                        q = "qnbinom({p}, size = {alpha}, prob = {(beta)} / (1 + {(beta)}){opts})", p_native = TRUE,
                        discrete = TRUE, support = c(0, Inf)),
    neg_binomial_2 = dist(c("mu", "phi"), r = "rnbinom({n}, size = {phi}, mu = {mu})",
                          d = "dnbinom({y}, size = {phi}, mu = {mu}, log = TRUE)",
                          p = "pnbinom({y}, size = {phi}, mu = {mu}{opts})",
                          q = "qnbinom({p}, size = {phi}, mu = {mu}{opts})", p_native = TRUE, discrete = TRUE,
                          support = c(0, Inf)),
    neg_binomial_2_log = dist(c("eta", "phi"), r = "rnbinom({n}, size = {phi}, mu = exp({eta}))",
                              d = "dnbinom({y}, size = {phi}, mu = exp({eta}), log = TRUE)",
                              p = "pnbinom({y}, size = {phi}, mu = exp({eta}){opts})",
                              q = "qnbinom({p}, size = {phi}, mu = exp({eta}){opts})", p_native = TRUE,
                              discrete = TRUE, support = c(0, Inf)),
    beta_neg_binomial = dist(c("r", "alpha", "beta"), r = "rnbinom({n}, size = {r}, prob = rbeta({n}, {alpha}, {beta}))",
                             d = "lchoose({y} + {(r)} - 1, {y}) + lbeta({(r)} + {(alpha)}, {y} + {(beta)}) - lbeta({alpha}, {beta})",
                             discrete = TRUE, support = c(0, Inf)),
    discrete_range = dist(c("lower", "upper"), r = "{(lower)} + floor(runif({n}) * ({(upper)} - {(lower)} + 1))",
                          d = "ifelse({y} >= {lower} & {y} <= {upper}, -log({(upper)} - {(lower)} + 1), -Inf)",
                          p = "pmin(1, pmax(0, (floor({y}) - {(lower)} + 1) / ({(upper)} - {(lower)} + 1)))",
                          q = "{(lower)} - 1 + ceiling({p} * ({(upper)} - {(lower)} + 1))", discrete = TRUE,
                          support = c("lower", "upper")),
    categorical = dist("theta", r = "rcategorical({n}, {theta})", d = "log({(theta)}[{y}])",
                       discrete = TRUE, helpers = "rcategorical"),
    categorical_logit = dist("beta", r = "rcategorical({n}, softmax({beta}))", d = "log_softmax({beta})[{y}]",
                             discrete = TRUE, helpers = c("rcategorical", "softmax", "log_softmax")),
    ordered_logistic = dist(c("eta", "c"), r = "rordered_logistic({eta}, {c})",
                            d = "dordered_logistic({y}, {eta}, {c})", discrete = TRUE,
                            helpers = c("rordered_logistic", "dordered_logistic")),
    ordered_probit = dist(c("eta", "c"), r = "rordered_probit({eta}, {c})",
                          d = "dordered_probit({y}, {eta}, {c})", discrete = TRUE,
                          helpers = c("rordered_probit", "dordered_probit")),
    # multivariate
    multinomial = dist("theta", r = "as.vector(rmultinom(1, size = sum({like}), prob = {theta}))",
                       d = "dmultinom({y}, prob = {theta}, log = TRUE)", discrete = TRUE, multivariate = TRUE,
                       note = "the number of trials is sum(outcome), so the outcome must be supplied in dat"),
    multinomial_logit = dist("gamma", r = "as.vector(rmultinom(1, size = sum({like}), prob = softmax({gamma})))",
                             d = "dmultinom({y}, prob = softmax({gamma}), log = TRUE)", discrete = TRUE,
                             multivariate = TRUE, helpers = "softmax",
                             note = "the number of trials is sum(outcome), so the outcome must be supplied in dat"),
    dirichlet = dist("alpha", r = "rdirichlet({like}, {alpha})", d = "ddirichlet({y}, {alpha})",
                     multivariate = TRUE, helpers = c("rdirichlet", "ddirichlet")),
    multi_normal = dist(c("mu", "Sigma"), r = "rmulti_normal({like}, {mu}, {Sigma})",
                        d = "dmulti_normal({y}, {mu}, {Sigma})", multivariate = TRUE,
                        helpers = c("rmulti_normal", "dmulti_normal")),
    multi_normal_cholesky = dist(c("mu", "L"), r = "rmulti_normal({like}, {mu}, tcrossprod({L}))",
                                 d = "dmulti_normal({y}, {mu}, tcrossprod({L}))", multivariate = TRUE,
                                 helpers = c("rmulti_normal", "dmulti_normal")),
    multi_normal_prec = dist(c("mu", "Omega"), r = "rmulti_normal({like}, {mu}, solve({Omega}))",
                             d = "dmulti_normal({y}, {mu}, solve({Omega}))", multivariate = TRUE,
                             helpers = c("rmulti_normal", "dmulti_normal")),
    multi_student_t = dist(c("nu", "mu", "Sigma"), r = "rmulti_student_t({like}, {nu}, {mu}, {Sigma})",
                           d = "dmulti_student_t({y}, {nu}, {mu}, {Sigma})", multivariate = TRUE,
                           helpers = c("rmulti_student_t", "dmulti_student_t")),
    multi_student_t_cholesky = dist(c("nu", "mu", "L"), r = "rmulti_student_t({like}, {nu}, {mu}, tcrossprod({L}))",
                                    d = "dmulti_student_t({y}, {nu}, {mu}, tcrossprod({L}))", multivariate = TRUE,
                                    helpers = c("rmulti_student_t", "dmulti_student_t")),
    lkj_corr = dist("eta", r = "rlkj_corr({like}, {eta})", d = "dlkj_corr({y}, {eta})", multivariate = TRUE,
                    helpers = c("rlkj_corr", "rlkj_corr_cholesky", "dlkj_corr")),
    lkj_corr_cholesky = dist("eta", r = "rlkj_corr_cholesky({like}, {eta})",
                             d = "dlkj_corr_cholesky({y}, {eta})", multivariate = TRUE,
                             helpers = c("rlkj_corr_cholesky", "dlkj_corr_cholesky", "dlkj_corr")),
    wishart = dist(c("nu", "Sigma"), r = "rWishart(1, {nu}, {Sigma})[, , 1]", d = "dwishart({y}, {nu}, {Sigma})",
                   multivariate = TRUE, helpers = "dwishart"),
    wishart_cholesky = dist(c("nu", "L_S"), r = "t(chol(rWishart(1, {nu}, tcrossprod({L_S}))[, , 1]))",
                            multivariate = TRUE),
    inv_wishart = dist(c("nu", "Sigma"), r = "solve(rWishart(1, {nu}, solve({Sigma}))[, , 1])",
                       d = "dinv_wishart({y}, {nu}, {Sigma})", multivariate = TRUE,
                       helpers = c("dinv_wishart", "dwishart")),
    inv_wishart_cholesky = dist(c("nu", "L_S"),
                                r = "t(chol(solve(rWishart(1, {nu}, solve(tcrossprod({L_S})))[, , 1])))",
                                multivariate = TRUE)
  )
}

# GLM forms: rewrite to their base distribution with a linear predictor
stan_glm_distributions <- c(normal_id_glm = "normal", bernoulli_logit_glm = "bernoulli_logit",
                            poisson_log_glm = "poisson_log", neg_binomial_2_log_glm = "neg_binomial_2_log")

# split a Stan density-function name into distribution and suffix
split_dist_name <- function(name) {
  m <- regmatches(name, regexec("^(.*)_(lpdf|lupdf|lpmf|lupmf|cdf|lcdf|lccdf|rng)$", name))[[1]]
  if (length(m) == 0) return(NULL)
  list(dist = m[2], suffix = m[3])
}

# substitute R code into a template
fill_template <- function(template, values) {
  out <- template
  for (nm in names(values)) {
    v <- values[[nm]]
    out <- gsub(paste0("{(", nm, ")}"), r_paren(v), out, fixed = TRUE)
    out <- gsub(paste0("{", nm, "}"), if (is.character(v)) v else v$code, out, fixed = TRUE)
  }
  left <- regmatches(out, gregexpr("\\{\\(?[A-Za-z_]+\\)?\\}", out))[[1]]
  if (length(left)) stop("Stan2R: internal error, unfilled template placeholder ", left[1], call. = FALSE)
  out
}

# values may be plain strings (treated as atomic) or translated expressions
r_paren <- function(v) {
  if (is.character(v)) return(v)
  if (v$prec >= PREC_ATOM) v$code else paste0("(", v$code, ")")
}
