// left-hand sides that are affine in one variable, as in the eGTEx spike-and-slab prior
data {
  int<lower=1> N;
  vector<lower=0>[N] s2;
  real<lower=0> nu;
}
parameters {
  real<lower=-1, upper=1> r;
  vector<lower=0, upper=1>[N] theta;
  real<lower=0> sigma;
}
model {
  (r + 1) / 2 ~ beta(2, 2);
  for (i in 1:N) {
    array[2] real lp;
    lp[1] = log(0.7) + normal_lpdf(theta[i] | 0.5, 0.2);
    lp[2] = log(0.3) + exponential_lpdf(1 - theta[i] | 8);
    target += log_sum_exp(lp);
  }
  sigma ~ lognormal(0, 0.3);
  nu * s2 / square(sigma) ~ chi_square(nu);
}
