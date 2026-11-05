dists <- Stan2R:::stan_distributions()

test_that("every distribution entry is well formed", {
  for (nm in names(dists)) {
    d <- dists[[nm]]
    expect_true(is.character(d$args), info = nm)
    for (w in c("r", "d", "p", "q")) {
      if (is.null(d[[w]])) next
      used <- regmatches(d[[w]], gregexpr("\\{\\(?[A-Za-z_]+\\)?\\}", d[[w]]))[[1]]
      used <- unique(gsub("[{}()]", "", used))
      expect_true(all(used %in% c(d$args, "n", "y", "p", "opts", "like")), info = paste(nm, w))
    }
    expect_true(all(d$helpers %in% names(Stan2R:::runtime_helpers)), info = nm)
  }
  expect_setequal(names(dist_test_args), names(dists)[!vapply(dists, `[[`, FALSE, "multivariate") &
                                                         !names(dists) %in% c("categorical", "categorical_logit", "ordered_logistic", "ordered_probit")])
})

test_that("continuous densities integrate to one and match their CDFs", {
  env <- dist_env()
  for (nm in names(dist_test_args)) {
    def <- dists[[nm]]
    if (def$discrete) next
    args <- dist_test_args[[nm]]
    d <- dist_fn(nm, "d", args, env)
    s <- dist_test_support(nm, def, args)
    total <- integrate(function(y) exp(d(y)), s[1], s[2], rel.tol = 1e-8, subdivisions = 1000)$value
    expect_equal(total, 1, tolerance = 1e-5, info = nm)
    if (!is.null(def$p)) {
      p <- dist_fn(nm, "p", args, env)
      lo <- if (is.finite(s[1])) s[1] else -Inf
      for (y in c(0.4, 1.3, 2.2)) {
        if (y <= s[1] || y >= s[2]) next
        expect_equal(p(y), integrate(function(t) exp(d(t)), lo, y, rel.tol = 1e-9)$value, tolerance = 1e-6, info = nm)
      }
    }
    if (!is.null(def$q)) {
      q <- dist_fn(nm, "q", args, env)
      p <- dist_fn(nm, "p", args, env)
      expect_equal(p(q(c(0.1, 0.5, 0.9))), c(0.1, 0.5, 0.9), tolerance = 1e-6, info = nm)
    }
  }
})

test_that("discrete masses sum to one and match their CDFs and quantiles", {
  env <- dist_env()
  for (nm in names(dist_test_args)) {
    def <- dists[[nm]]
    if (!def$discrete) next
    args <- dist_test_args[[nm]]
    d <- dist_fn(nm, "d", args, env)
    s <- dist_test_support(nm, def, args)
    ys <- s[1]:min(s[2], 200)
    pm <- exp(d(ys))
    expect_equal(sum(pm), 1, tolerance = 1e-8, info = nm)
    if (!is.null(def$p)) expect_equal(dist_fn(nm, "p", args, env)(ys), cumsum(pm), tolerance = 1e-8, info = nm)
    if (!is.null(def$q)) {
      q <- dist_fn(nm, "q", args, env)
      cdf <- cumsum(pm)
      probe <- c(0.07, 0.5, 0.93)
      expect_equal(q(probe), vapply(probe, function(u) ys[which(cdf >= u)[1]], 0), info = nm)
    }
  }
})

test_that("random draws follow the density (KS / chi-square)", {
  env <- dist_env()
  set.seed(123)
  for (nm in names(dist_test_args)) {
    def <- dists[[nm]]
    args <- dist_test_args[[nm]]
    r <- dist_fn(nm, "r", args, env)
    x <- r(20000)
    d <- dist_fn(nm, "d", args, env)
    if (def$discrete) {
      s <- dist_test_support(nm, def, args)
      ys <- s[1]:min(s[2], 200)
      expected <- exp(d(ys)) * length(x)
      keep <- expected > 5
      observed <- vapply(ys[keep], function(k) sum(x == k), 0)
      stat <- sum((observed - expected[keep])^2 / expected[keep])
      expect_lt(pchisq(stat, sum(keep) - 1, lower.tail = FALSE), 1, label = nm)
      expect_gt(pchisq(stat, sum(keep) - 1, lower.tail = FALSE), 1e-4, label = nm)
    } else if (!is.null(def$p)) {
      p <- dist_fn(nm, "p", args, env)
      expect_gt(suppressWarnings(ks.test(x, p)$p.value), 1e-4, label = nm)
    } else {
      s <- dist_test_support(nm, def, args)
      m <- integrate(function(y) y * exp(d(y)), s[1], s[2])$value
      expect_equal(mean(x), m, tolerance = 0.05, label = nm)
    }
  }
})

test_that("multivariate helpers produce the right shapes and moments", {
  env <- dist_env()
  set.seed(7)
  Sigma <- matrix(c(2, 0.6, 0.6, 1), 2)
  draws <- env$rmulti_normal(matrix(NA, 5000, 2), c(1, -1), Sigma)
  expect_equal(dim(draws), c(5000, 2))
  expect_equal(colMeans(draws), c(1, -1), tolerance = 0.05)
  expect_equal(cov(draws), Sigma, tolerance = 0.06)
  y <- c(0.3, -0.2)
  expect_equal(env$dmulti_normal(y, c(1, -1), Sigma),
               -log(2 * pi) - 0.5 * log(det(Sigma)) - 0.5 * sum((y - c(1, -1)) * solve(Sigma, y - c(1, -1))))
  L <- env$rlkj_corr_cholesky(diag(4), 2)
  expect_equal(diag(tcrossprod(L)), rep(1, 4))
  expect_true(all(L[upper.tri(L)] == 0))
  expect_equal(sum(env$rdirichlet(rep(NA, 3), c(1, 2, 3))), 1)
  # LKJ(eta) marginal of one correlation is 2 * Beta(eta - 1 + K/2, same) - 1
  r12 <- replicate(4000, tcrossprod(env$rlkj_corr_cholesky(diag(3), 2))[1, 2])
  expect_gt(ks.test((r12 + 1) / 2, function(q) pbeta(q, 2.5, 2.5))$p.value, 1e-4)
})
