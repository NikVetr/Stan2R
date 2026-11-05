data {
  int<lower=0> N;
  int<lower=1> K;
  matrix[N, K] X;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
  vector[N] y_raw;  // outcome moved from data, non-centred: y = mu + sigma * y_raw
}
transformed parameters {
  vector[N] mu = alpha + X * beta;
  vector[N] y = mu + sigma * y_raw;
}
model {
  alpha ~ normal(0, 2);
  beta ~ normal(0, 1);
  sigma ~ student_t(3, 0, 2.5);
  y_raw ~ std_normal();
}
