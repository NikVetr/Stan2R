// Poisson regression with an exposure offset and a sum-to-zero group effect
data {
  int<lower=1> N;
  int<lower=2> J;
  array[N] int<lower=1, upper=J> g;
  vector[N] log_exposure;
  array[N] int<lower=0> y;
}
parameters {
  real alpha;
  sum_to_zero_vector[J] u;
}
model {
  alpha ~ normal(1, 0.5);
  u ~ normal(0, 0.7);
  y ~ poisson_log(alpha + u[g] + log_exposure);
}
