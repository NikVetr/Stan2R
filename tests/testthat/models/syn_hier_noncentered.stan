// non-centred hierarchical model with a transformed-parameter scale
data {
  int<lower=1> J;
  int<lower=1> N;
  array[N] int<lower=1, upper=J> g;
  vector[N] y;
}
parameters {
  real mu;
  real log_tau;
  vector[J] z;
  real<lower=0> sigma;
}
transformed parameters {
  real<lower=0> tau = exp(log_tau);
  vector[J] theta = mu + tau * z;
}
model {
  mu ~ normal(0, 5);
  log_tau ~ normal(-1, 0.5);
  z ~ std_normal();
  sigma ~ lognormal(0, 0.5);
  for (n in 1:N) {
    y[n] ~ normal(theta[g[n]], sigma);
  }
}
