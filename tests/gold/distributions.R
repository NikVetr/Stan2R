# Gold standard for the distribution table: Stan's own _lpdf/_lpmf, _lcdf and
# _rng (via fixed_param) against Stan2R's d, p and r templates.
#
#   Rscript tests/gold/distributions.R [out_dir]     (compiles one model; 4 cores max)

out_dir <- commandArgs(TRUE)[1]
if (is.na(out_dir)) out_dir <- "tmp/output/gold"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-distributions.R")
library(cmdstanr)

dists <- Stan2R:::stan_distributions()
math <- system2(Stan2R:::stanc_path(), "--dump-stan-math-distributions", stdout = TRUE)
suffixes <- setNames(strsplit(sub("^[^:]+: ", "", math), ", "), sub(":.*$", "", math))

points <- function(nm, def, args) {
  s <- dist_test_support(nm, def, args)
  if (def$discrete) return(unique(pmin(pmax(c(s[1], s[1] + 1, s[1] + 3), s[1]), s[2])))
  lo <- if (is.finite(s[1])) s[1] else -2
  hi <- if (is.finite(s[2])) s[2] else lo + 6
  lo + (hi - lo) * c(0.13, 0.41, 0.77)
}
num <- function(x) format(x, digits = 17)
# Stan arguments: integers where Stan requires them, otherwise real literals
int_args <- list(binomial = 1, binomial_logit = 1, beta_binomial = 1, hypergeometric = 1:3, discrete_range = 1:2)
arg_literals <- function(nm, args) {
  vapply(seq_along(args), function(i) {
    if (i %in% int_args[[nm]]) as.character(as.integer(args[[i]]))
    else { v <- num(args[[i]]); if (grepl("^-?[0-9]+$", v)) paste0(v, ".0") else v }
  }, "")
}

# which CDF form stanc accepts for each distribution (the dump's suffix list is not exact)
stanc_accepts <- function(line) {
  f <- tempfile(fileext = ".stan")
  writeLines(c("generated quantities {", line, "}"), f)
  system2(Stan2R:::stanc_path(), c(shQuote(f), "-o", tempfile()), stdout = FALSE, stderr = FALSE) == 0
}
cdf_fun <- list()
for (nm in names(dist_test_args)) {
  def <- dists[[nm]]
  args <- dist_test_args[[nm]]
  arg_text <- paste(arg_literals(nm, args), collapse = ", ")
  y <- if (def$discrete) "1" else "0.5"
  sep <- if (nzchar(arg_text)) " | " else ""
  cdf_fun[[nm]] <- if (stanc_accepts(sprintf("  real c = %s_cdf(%s%s%s);", nm, y, sep, arg_text))) "cdf"
                   else if (stanc_accepts(sprintf("  real c = %s_lcdf(%s%s%s);", nm, y, sep, arg_text))) "lcdf" else ""
}

gq <- character(0)
checks <- list()
for (nm in names(dist_test_args)) {
  def <- dists[[nm]]
  args <- dist_test_args[[nm]]
  sx <- suffixes[[nm]]
  if (is.null(sx)) next
  arg_text <- paste(arg_literals(nm, args), collapse = ", ")
  lp <- if (def$discrete) "lpmf" else "lpdf"
  ys <- points(nm, def, args)
  for (k in seq_along(ys)) {
    y <- if (def$discrete) as.character(as.integer(ys[k])) else num(ys[k])
    sep <- if (nzchar(arg_text)) " | " else ""
    gq <- c(gq, sprintf("  real lp_%s_%d = %s_%s(%s%s%s);", nm, k, nm, lp, y, sep, arg_text))
    if (cdf_fun[[nm]] != "") gq <- c(gq, sprintf("  real cdf_%s_%d = %s(%s_%s(%s%s%s));", nm, k,
                                                  if (cdf_fun[[nm]] == "lcdf") "exp" else "", nm, cdf_fun[[nm]], y, sep, arg_text))
  }
  if ("rng" %in% sx) gq <- c(gq, sprintf("  %s rng_%s = %s_rng(%s);", if (def$discrete) "int" else "real", nm, nm, arg_text))
  checks[[nm]] <- list(ys = ys, cdf = cdf_fun[[nm]] != "", rng = "rng" %in% sx)
}
model_file <- file.path(out_dir, "distribution_check.stan")
writeLines(c("generated quantities {", gq, "}"), model_file)

mod <- cmdstan_model(model_file, quiet = TRUE)
fit <- mod$sample(fixed_param = TRUE, chains = 1, iter_sampling = 5000, seed = 20261104, refresh = 0, threads_per_chain = 1,
                  show_messages = FALSE, output_dir = out_dir)
draws <- posterior::as_draws_df(fit$draws())

env <- dist_env()
rows <- list()
set.seed(20261104)
for (nm in names(checks)) {
  ck <- checks[[nm]]
  args <- dist_test_args[[nm]]
  d <- dist_fn(nm, "d", args, env)
  for (k in seq_along(ck$ys)) {
    stan_lp <- draws[[sprintf("lp_%s_%d", nm, k)]][1]
    rows[[length(rows) + 1]] <- data.frame(dist = nm, check = sprintf("log density at %s", signif(ck$ys[k], 4)),
                                           stan = stan_lp, stan2r = d(ck$ys[k]), stat = abs(stan_lp - d(ck$ys[k])),
                                           pass = isTRUE(all.equal(stan_lp, d(ck$ys[k]), tolerance = 1e-7)))
    if (ck$cdf && !is.null(dists[[nm]]$p)) {
      stan_c <- draws[[sprintf("cdf_%s_%d", nm, k)]][1]
      p <- dist_fn(nm, "p", args, env)(ck$ys[k])
      rows[[length(rows) + 1]] <- data.frame(dist = nm, check = sprintf("CDF at %s", signif(ck$ys[k], 4)),
                                             stan = stan_c, stan2r = p, stat = abs(stan_c - p),
                                             pass = abs(stan_c - p) < 1e-7)
    }
  }
  if (ck$rng) {
    s <- draws[[paste0("rng_", nm)]]
    n_failed <- sum(!is.finite(s))
    s <- s[is.finite(s)]
    r <- dist_fn(nm, "r", args, env)(length(s))
    if (dists[[nm]]$discrete) {
      lv <- sort(unique(c(s, r)))
      tab <- rbind(table(factor(s, lv)), table(factor(r, lv)))
      tab <- tab[, colSums(tab) >= 10, drop = FALSE]
      pval <- if (ncol(tab) > 1) suppressWarnings(chisq.test(tab)$p.value) else 1
      label <- "draws vs Stan rng (chi-square p)"
    } else {
      pval <- suppressWarnings(ks.test(s, r)$p.value)
      label <- "draws vs Stan rng (KS p)"
    }
    if (n_failed) label <- sprintf("%s; %d Stan draws failed", label, n_failed)
    rows[[length(rows) + 1]] <- data.frame(dist = nm, check = label, stan = mean(s), stan2r = mean(r), stat = pval,
                                           pass = pval > 1e-3)
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "distributions.csv"), row.names = FALSE)
cat(sprintf("%d checks over %d distributions; %d failed\n", nrow(res), length(unique(res$dist)), sum(!res$pass)))
print(res[!res$pass, ], row.names = FALSE)
