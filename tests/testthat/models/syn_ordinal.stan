// ordinal regression with ordered cutpoints
data {
  int<lower=1> N;
  int<lower=3> K;
  vector[N] x;
  array[N] int<lower=1, upper=K> y;
}
parameters {
  real beta;
  ordered[K - 1] c;
}
model {
  beta ~ normal(0, 1);
  c ~ normal(0, 2);
  y ~ ordered_logistic(x * beta, c);
}
