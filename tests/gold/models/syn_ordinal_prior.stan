data {
  int<lower=1> N;
  int<lower=3> K;
  vector[N] x;
}
parameters {
  real beta;
  ordered[K - 1] c;
}
model {
  beta ~ normal(0, 1);
  c ~ normal(0, 2);
}
generated quantities {
  array[N] int y_sim;
  for (n in 1:N) y_sim[n] = ordered_logistic_rng(x[n] * beta, c);
}
