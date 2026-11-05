data {
  int<lower=1> N;
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
}
generated quantities {
  vector[N] y_sim;
  for (n in 1:N) {
    y_sim[n] = bernoulli_rng(lambda) ? normal_rng(mu[1], sigma) : normal_rng(mu[2], sigma);
  }
}
