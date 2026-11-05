// mixture prior on a bounded parameter: each component lives on [0, 1]
data {
  int<lower=1> N;
  array[N] int<lower=0> trials;
  array[N] int<lower=0> successes;
  real<lower=0> s;
  real<lower=0> rate;
}
parameters {
  vector<lower=0, upper=1>[N] theta;
}
model {
  for (i in 1:N) {
    target += log_sum_exp(log(0.6) + normal_lpdf(theta[i] | 0.5, s), log(0.4) + exponential_lpdf(theta[i] | rate));
  }
  successes ~ binomial(trials, theta);
}
