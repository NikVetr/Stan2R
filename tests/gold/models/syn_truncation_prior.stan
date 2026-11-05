data {
  int<lower=1> N;
  real<lower=0> U;
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
  vector<lower=0, upper=U>[N] y;  // outcome moved from data, on its truncation interval
}
model {
  mu ~ double_exponential(1, 2);
  sigma ~ gamma(4, 4);
  for (n in 1:N) {
    y[n] ~ normal(mu, sigma) T[0, U];
  }
}
