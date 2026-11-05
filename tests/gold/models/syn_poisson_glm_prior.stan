data {
  int<lower=1> N;
  int<lower=2> J;
  array[N] int<lower=1, upper=J> g;
  vector[N] log_exposure;
}
parameters {
  real alpha;
  sum_to_zero_vector[J] u;
}
model {
  alpha ~ normal(1, 0.5);
  u ~ normal(0, 0.7);
}
generated quantities {
  array[N] int y_sim = poisson_log_rng(alpha + u[g] + log_exposure);
}
