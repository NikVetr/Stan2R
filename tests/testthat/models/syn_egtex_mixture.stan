// three-component allele-count mixture written with an lp accumulator (as in the eGTEx ASE models)
data {
  int<lower=1> n;
  array[n] int<lower=1> total;
  array[n] int<lower=0> count;
}
parameters {
  simplex[3] w;
  real<lower=0> conc;
  real<lower=0, upper=1> loc;
}
transformed parameters {
  real a = loc * conc;
  real b = (1 - loc) * conc;
}
model {
  w ~ dirichlet(rep_vector(2, 3));
  conc ~ lognormal(2, 0.5);
  loc ~ beta(3, 3);
  for (i in 1:n) {
    array[3] real lp;
    lp[1] = log(w[1]) + binomial_lpmf(count[i] | total[i], 0.5);
    lp[2] = log(w[2]) + beta_binomial_lpmf(count[i] | total[i], a, b);
    lp[3] = log(w[3]) + beta_binomial_lpmf(count[i] | total[i], b, a);
    target += log_sum_exp(lp);
  }
}
