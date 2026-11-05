# Runtime helpers for generated R code. Each entry is the source text of one
# base-R function; parse_Stan() copies the ones a model uses into the top of
# the generated script so the script runs on its own (no Stan2R needed).
# `needs` lists helpers an entry calls.

runtime_helpers <- list(
  stan_range = list(code = '
# Stan range a:b, which is empty when b < a
stan_range <- function(a, b) if (b >= a) a:b else integer(0)'),

  stan_modulus = list(code = '
# Stan (C++) integer modulus: the result takes the sign of the dividend
stan_modulus <- function(a, b) a - b * trunc(a / b)'),

  stan_round = list(code = '
# Stan rounds halves away from zero (R rounds them to even)
stan_round <- function(x) sign(x) * floor(abs(x) + 0.5)'),

  drop_dims = list(code = '
# drop the array dimensions that a single index removed in Stan
drop_dims <- function(x, dims) {
  d <- dim(x)[-dims]
  if (length(d) <= 1) as.vector(x) else array(x, dim = d)
}'),

  stack_first_dim = list(code = '
# Stan array expression {a, b, ...} of containers: stack along a new first dimension
stack_first_dim <- function(elements) {
  inner <- dim(elements[[1]])
  if (is.null(inner)) inner <- length(elements[[1]])
  out <- array(NA_real_, dim = c(length(elements), inner))
  for (i in seq_along(elements)) out[i + length(elements) * (seq_len(prod(inner)) - 1)] <- elements[[i]]
  out
}'),

  log_sum_exp = list(code = '
log_sum_exp <- function(...) {
  x <- c(...)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}'),

  log_mix = list(code = '
# log_mix(theta, lp1, lp2) for a scalar mixing proportion, or
# log_mix(theta_vector, lp_vector) for K components
log_mix <- function(theta, lp1, lp2 = NULL) {
  if (is.null(lp2)) return(log_sum_exp(log(theta) + lp1))
  log_sum_exp(log(theta) + lp1, log1p(-theta) + lp2)
}', needs = "log_sum_exp"),

  softmax = list(code = '
softmax <- function(x) {
  e <- exp(x - max(x))
  e / sum(e)
}'),

  log_softmax = list(code = '
log_softmax <- function(x) x - log_sum_exp(x)', needs = "log_sum_exp"),

  add_log_density = list(code = '
# accumulate pointwise log densities; NA marks elements without a term yet
add_log_density <- function(current, term) ifelse(is.na(current), 0, current) + term'),

  plaplace = list(code = '
plaplace <- function(q, mu, sigma) ifelse(q < mu, 0.5 * exp((q - mu) / sigma), 1 - 0.5 * exp(-(q - mu) / sigma))'),

  qlaplace = list(code = '
qlaplace <- function(p, mu, sigma) ifelse(p < 0.5, mu + sigma * log(2 * p), mu - sigma * log(2 - 2 * p))'),

  dexp_mod_normal = list(code = '
dexp_mod_normal <- function(y, mu, sigma, lambda) {
  log(lambda) + lambda * (mu + 0.5 * lambda * sigma^2 - y) +
    pnorm((y - mu) / sigma - lambda * sigma, log.p = TRUE)
}'),

  pexp_mod_normal = list(code = '
pexp_mod_normal <- function(y, mu, sigma, lambda) {
  pnorm(y, mu, sigma) - exp(-lambda * (y - mu) + 0.5 * (lambda * sigma)^2) *
    pnorm((y - mu) / sigma - lambda * sigma)
}'),

  rskew_normal = list(code = '
rskew_normal <- function(n, xi, omega, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  xi + omega * (delta * abs(rnorm(n)) + sqrt(1 - delta^2) * rnorm(n))
}'),

  pbeta_binomial = list(code = '
pbeta_binomial <- function(q, N, alpha, beta) {
  mapply(function(q, N, alpha, beta) {
    if (q < 0) return(0)
    k <- 0:min(floor(q), N)
    sum(exp(lchoose(N, k) + lbeta(k + alpha, N - k + beta) - lbeta(alpha, beta)))
  }, q, N, alpha, beta)
}'),

  qbeta_binomial = list(code = '
qbeta_binomial <- function(p, N, alpha, beta) {
  mapply(function(p, N, alpha, beta) {
    k <- 0:N
    cdf <- cumsum(exp(lchoose(N, k) + lbeta(k + alpha, N - k + beta) - lbeta(alpha, beta)))
    k[which(cdf >= p - 1e-12)[1]]
  }, p, N, alpha, beta)
}'),

  rcategorical = list(code = '
rcategorical <- function(n, theta) sample.int(length(theta), n, replace = TRUE, prob = theta)'),

  rordered_logistic = list(code = '
# P(y = k) = inv_logit(eta - c[k - 1]) - inv_logit(eta - c[k]), with c[0] = -Inf, c[K] = Inf
rordered_logistic <- function(eta, c) {
  vapply(eta, function(e) {
    p <- -diff(c(1, plogis(e - c), 0))
    sample.int(length(p), 1, prob = p)
  }, numeric(1))
}'),

  dordered_logistic = list(code = '
dordered_logistic <- function(y, eta, c) {
  eta <- rep_len(eta, length(y))
  vapply(seq_along(y), function(i) log(-diff(c(1, plogis(eta[i] - c), 0)))[y[i]], numeric(1))
}'),

  rordered_probit = list(code = '
rordered_probit <- function(eta, c) {
  vapply(eta, function(e) {
    p <- -diff(c(1, pnorm(e - c), 0))
    sample.int(length(p), 1, prob = p)
  }, numeric(1))
}'),

  dordered_probit = list(code = '
dordered_probit <- function(y, eta, c) {
  eta <- rep_len(eta, length(y))
  vapply(seq_along(y), function(i) log(-diff(c(1, pnorm(eta[i] - c), 0)))[y[i]], numeric(1))
}'),

  rdirichlet = list(code = '
# fill `like` (a simplex, or an array of simplexes stored one per row)
rdirichlet <- function(like, alpha) {
  one <- function(a) {
    g <- rgamma(length(a), a)
    g / sum(g)
  }
  if (is.matrix(like)) {
    alpha <- matrix(alpha, nrow(like), ncol(like), byrow = !is.matrix(alpha))
    t(apply(alpha, 1, one))
  } else {
    one(rep_len(alpha, length(like)))
  }
}'),

  ddirichlet = list(code = '
ddirichlet <- function(y, alpha) {
  one <- function(y, a) lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * log(y))
  if (is.matrix(y)) {
    alpha <- matrix(alpha, nrow(y), ncol(y), byrow = !is.matrix(alpha))
    vapply(seq_len(nrow(y)), function(i) one(y[i, ], alpha[i, ]), numeric(1))
  } else {
    one(y, alpha)
  }
}'),

  rmulti_normal = list(code = '
# fill `like` (a vector, or an array of vectors stored one per row) with draws
rmulti_normal <- function(like, mu, Sigma) {
  R <- chol(Sigma)
  if (is.matrix(like)) {
    z <- matrix(rnorm(length(like)), nrow(like), ncol(like)) %*% R
    if (is.matrix(mu)) z + mu else sweep(z, 2, mu, "+")
  } else {
    as.vector(mu + rnorm(length(like)) %*% R)
  }
}'),

  dmulti_normal = list(code = '
# log density per vector (per row when y stores an array of vectors)
dmulti_normal <- function(y, mu, Sigma) {
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1)
  mu <- if (is.matrix(mu)) mu else matrix(mu, nrow(y), ncol(y), byrow = TRUE)
  R <- chol(Sigma)
  z <- forwardsolve(t(R), t(y - mu))
  -0.5 * ncol(y) * log(2 * pi) - sum(log(diag(R))) - 0.5 * colSums(z^2)
}'),

  rmulti_student_t = list(code = '
rmulti_student_t <- function(like, nu, mu, Sigma) {
  R <- chol(Sigma)
  if (is.matrix(like)) {
    w <- sqrt(nu / rchisq(nrow(like), nu))
    z <- (matrix(rnorm(length(like)), nrow(like), ncol(like)) %*% R) * w
    if (is.matrix(mu)) z + mu else sweep(z, 2, mu, "+")
  } else {
    as.vector(mu + sqrt(nu / rchisq(1, nu)) * rnorm(length(like)) %*% R)
  }
}'),

  dmulti_student_t = list(code = '
dmulti_student_t <- function(y, nu, mu, Sigma) {
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1)
  mu <- if (is.matrix(mu)) mu else matrix(mu, nrow(y), ncol(y), byrow = TRUE)
  K <- ncol(y)
  R <- chol(Sigma)
  z <- forwardsolve(t(R), t(y - mu))
  lgamma((nu + K) / 2) - lgamma(nu / 2) - 0.5 * K * log(nu * pi) - sum(log(diag(R))) -
    0.5 * (nu + K) * log1p(colSums(z^2) / nu)
}'),

  rlkj_corr_cholesky = list(code = '
# onion method (Lewandowski, Kurowicka and Joe 2009); returns the lower Cholesky factor
rlkj_corr_cholesky <- function(like, eta) {
  K <- nrow(like)
  if (K == 1) return(matrix(1, 1, 1))
  alpha <- eta + (K - 2) / 2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) {
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)
      z <- rnorm(m)
      R[1:m, m + 1] <- sqrt(y) * z / sqrt(sum(z^2))
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  }
  t(R)
}'),

  rlkj_corr = list(code = '
rlkj_corr <- function(like, eta) tcrossprod(rlkj_corr_cholesky(like, eta))', needs = "rlkj_corr_cholesky"),

  dlkj_corr = list(code = '
# LKJ log density including its normalising constant
dlkj_corr <- function(y, eta) {
  K <- nrow(y)
  log_c <- 0
  for (i in seq_len(K - 1)) {
    b <- eta + (K - 1 - i) / 2
    log_c <- log_c + (2 * eta - 2 + K - i) * (K - i) * log(2) +
      (K - i) * lbeta(b, b)
  }
  (eta - 1) * as.numeric(determinant(y)$modulus) - log_c
}'),

  dlkj_corr_cholesky = list(code = '
dlkj_corr_cholesky <- function(y, eta) {
  K <- nrow(y)
  dlkj_corr(tcrossprod(y), eta) + sum((K - seq_len(K)) * log(diag(y)))
}', needs = "dlkj_corr"),

  dwishart = list(code = '
dwishart <- function(W, nu, Sigma) {
  K <- nrow(W)
  ld <- function(m) as.numeric(determinant(m)$modulus)
  0.5 * (nu - K - 1) * ld(W) - 0.5 * sum(diag(solve(Sigma, W))) - 0.5 * nu * K * log(2) -
    0.5 * nu * ld(Sigma) - (0.25 * K * (K - 1) * log(pi) + sum(lgamma((nu + 1 - seq_len(K)) / 2)))
}'),

  dinv_wishart = list(code = '
dinv_wishart <- function(W, nu, Sigma) {
  K <- nrow(W)
  dwishart(solve(W), nu, solve(Sigma)) - (K + 1) * as.numeric(determinant(W)$modulus)
}', needs = "dwishart"),

  rtruncated_reject = list(code = '
# draws from a distribution restricted to [lower, upper] by rejection
rtruncated_reject <- function(n, draw, lower, upper, max_tries = 1e5) {
  out <- draw(n)
  bad <- out < lower | out > upper
  tries <- 0
  while (any(bad)) {
    tries <- tries + 1
    if (tries > max_tries) stop("rtruncated_reject: no draw inside the bounds after ", max_tries, " tries")
    out[bad] <- draw(n)[bad]
    bad <- out < lower | out > upper
  }
  out
}'),

  qr_thin_Q = list(code = '
# thin QR with the sign convention of Stan (R has a non-negative diagonal)
qr_thin_Q <- function(A) {
  d <- qr(A)
  qr.Q(d) %*% diag(sign(diag(qr.R(d))) + (diag(qr.R(d)) == 0), nrow = min(dim(A)))
}'),

  qr_thin_R = list(code = '
qr_thin_R <- function(A) {
  d <- qr(A)
  diag(sign(diag(qr.R(d))) + (diag(qr.R(d)) == 0), nrow = min(dim(A))) %*% qr.R(d)
}'),

  qr_Q = list(code = '
qr_Q <- function(A) {
  d <- qr(A)
  Q <- qr.Q(d, complete = TRUE)
  s <- sign(diag(qr.R(d, complete = TRUE)))
  s[s == 0 | is.na(s)] <- 1
  Q[, seq_along(s)] <- Q[, seq_along(s)] %*% diag(s, nrow = length(s))
  Q
}'),

  qr_R = list(code = '
qr_R <- function(A) {
  d <- qr(A)
  R <- qr.R(d, complete = TRUE)
  s <- sign(diag(R))
  s[s == 0 | is.na(s)] <- 1
  R[seq_along(s), ] <- diag(s, nrow = length(s)) %*% R[seq_along(s), , drop = FALSE]
  R
}'),

  matrix_power = list(code = '
matrix_power <- function(m, k) {
  out <- diag(nrow(m))
  for (i in seq_len(k)) out <- out %*% m
  out
}')
)

# resolve helper dependencies and return their source in a stable order
helper_source <- function(names) {
  todo <- unique(names)
  seen <- character(0)
  while (length(todo)) {
    h <- todo[1]
    todo <- todo[-1]
    if (h %in% seen) next
    if (is.null(runtime_helpers[[h]])) stop("Stan2R: internal error, unknown runtime helper ", h, call. = FALSE)
    seen <- c(seen, h)
    todo <- c(todo, runtime_helpers[[h]]$needs)
  }
  ordered <- names(runtime_helpers)[names(runtime_helpers) %in% seen]
  vapply(ordered, function(h) sub("^\n", "", runtime_helpers[[h]]$code), "")
}
