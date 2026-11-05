data {
  int<lower=1> n;
  array[n] int<lower=1> total;
}
parameters {
  vector<lower=0, upper=1>[n] theta;
  real<lower=0> lambda;
}
model {
  lambda ~ lognormal(1, 0.3);
  real tail_mass = exponential_lcdf(1 | lambda);
  for (i in 1:n) {
    array[3] real lp_theta;
    lp_theta[1] = log(0.5) + normal_lpdf(theta[i] | 0.5, 0.15);
    lp_theta[2] = log(0.25) + exponential_lpdf(theta[i] | lambda) - tail_mass;
    lp_theta[3] = log(0.25) + exponential_lpdf(1 - theta[i] | lambda) - tail_mass;
    target += log_sum_exp(lp_theta);
  }
}
generated quantities {
  array[n] int count_sim;
  for (i in 1:n) count_sim[i] = bernoulli_rng(0.4) ? binomial_rng(total[i], 0.5) : binomial_rng(total[i], theta[i]);
}
