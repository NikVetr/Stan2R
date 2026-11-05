// user-defined functions: a helper, a density with its own _rng, and a density without one
functions {
  real scaled_logit(real x, real s) {
    return inv_logit(x / s);
  }
  real shifted_exp_lpdf(real y, real shift, real rate) {
    return exponential_lpdf(y - shift | rate);
  }
  real shifted_exp_rng(real shift, real rate) {
    return shift + exponential_rng(rate);
  }
  real two_normals_lpdf(real y, real w, real m) {
    return log_sum_exp(log(w) + normal_lpdf(y | -m, 1), log1m(w) + normal_lpdf(y | m, 1));
  }
}
data {
  int<lower=1> N;
  vector[N] y;
  array[N] int<lower=0, upper=1> b;
}
parameters {
  real a;
  real<lower=0> m;
}
model {
  a ~ normal(0, 1);
  m ~ shifted_exp(0.5, 2);
  for (n in 1:N) {
    y[n] ~ two_normals(0.3, m);
    b[n] ~ bernoulli(scaled_logit(a, 2));
  }
}
