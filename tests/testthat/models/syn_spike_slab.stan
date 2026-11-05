// per-observation mixture likelihood that uses a parameter whose mixture prior
// is written later in the same loop (the layout of the eGTEx spike-and-slab models)
data {
  int<lower=1> n;
  array[n] int<lower=1> total;
  array[n] int<lower=0> count;
}
parameters {
  vector<lower=0, upper=1>[n] theta;
  real<lower=0> lambda;
}
model {
  lambda ~ lognormal(1, 0.3);
  real tail_mass = exponential_lcdf(1 | lambda);
  for (i in 1:n) {
    // likelihood mixture
    array[2] real lp_count;
    lp_count[1] = log(0.4) + binomial_lpmf(count[i] | total[i], 0.5);
    lp_count[2] = log(0.6) + binomial_lpmf(count[i] | total[i], theta[i]);
    target += log_sum_exp(lp_count);
    // theta prior: slab around 0.5 or a spike near either end
    array[3] real lp_theta;
    lp_theta[1] = log(0.5) + normal_lpdf(theta[i] | 0.5, 0.15);
    lp_theta[2] = log(0.25) + exponential_lpdf(theta[i] | lambda) - tail_mass;
    lp_theta[3] = log(0.25) + exponential_lpdf(1 - theta[i] | lambda) - tail_mass;
    target += log_sum_exp(lp_theta);
  }
}
