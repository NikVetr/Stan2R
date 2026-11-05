data {
  int<lower=1> N;
  int<lower=2> K;
}
parameters {
  vector[K] mu;
  vector<lower=0>[K] tau;
  cholesky_factor_corr[K] L;
}
transformed parameters {
  matrix[K, K] L_Sigma = diag_pre_multiply(tau, L);
}
model {
  mu ~ normal(0, 1);
  tau ~ cauchy(0, 1);
  L ~ lkj_corr_cholesky(2);
}
generated quantities {
  matrix[K, K] Omega = multiply_lower_tri_self_transpose(L);
  array[N] vector[K] y_sim;
  for (n in 1:N) y_sim[n] = multi_normal_cholesky_rng(mu, L_Sigma);
}
