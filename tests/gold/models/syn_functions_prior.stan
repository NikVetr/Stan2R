functions {
  real scaled_logit(real x, real s) {
    return inv_logit(x / s);
  }
  real shifted_exp_lpdf(real y, real shift, real rate) {
    return exponential_lpdf(y - shift | rate);
  }
  real two_normals_lpdf(real y, real w, real m) {
    return log_sum_exp(log(w) + normal_lpdf(y | -m, 1), log1m(w) + normal_lpdf(y | m, 1));
  }
}
data {
  int<lower=1> N;
}
parameters {
  real a;
  real<lower=0.5> m;
  vector[N] y;  // outcome moved from data
}
model {
  a ~ normal(0, 1);
  m ~ shifted_exp(0.5, 2);
  for (n in 1:N) {
    y[n] ~ two_normals(0.3, m);
  }
}
generated quantities {
  array[N] int b_sim;
  for (n in 1:N) b_sim[n] = bernoulli_rng(scaled_logit(a, 2));
}
