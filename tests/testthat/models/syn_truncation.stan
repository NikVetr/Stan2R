// explicit truncation of a likelihood, a bounded Laplace prior, and while/if control flow
data {
  int<lower=1> N;
  vector<lower=0>[N] y;
  real<lower=0> U;
}
transformed data {
  int n_half = N %/% 2;
  real total = 0;
  int i = 1;
  while (i <= N) {
    if (i > n_half) {
      total += 1;
    } else {
      total += 0.5;
    }
    i += 1;
  }
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
}
model {
  mu ~ double_exponential(1, 2);
  sigma ~ gamma(4, 4);
  for (n in 1:N) {
    y[n] ~ normal(mu, sigma) T[0, U];
  }
}
