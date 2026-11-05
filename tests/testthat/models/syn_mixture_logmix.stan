// two-component normal mixture written with log_mix and ordered means
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  real<lower=0, upper=1> lambda;
  ordered[2] mu;
  real<lower=0> sigma;
}
model {
  lambda ~ beta(2, 2);
  mu ~ normal(0, 3);
  sigma ~ exponential(1);
  for (n in 1:N) {
    target += log_mix(lambda, normal_lpdf(y[n] | mu[1], sigma), normal_lpdf(y[n] | mu[2], sigma));
  }
}
