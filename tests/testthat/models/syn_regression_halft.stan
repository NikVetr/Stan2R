// linear regression with a matrix of predictors and a half-t scale prior
data {
  int<lower=0> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] mu = alpha + X * beta;
}
model {
  alpha ~ normal(0, 2);
  beta ~ normal(0, 1);
  sigma ~ student_t(3, 0, 2.5);
  y ~ normal(mu, sigma);
}
generated quantities {
  array[N] real y_rep = normal_rng(mu, sigma);
}
