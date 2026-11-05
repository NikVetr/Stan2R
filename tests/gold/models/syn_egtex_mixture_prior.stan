data {
  int<lower=1> n;
  array[n] int<lower=1> total;
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
}
generated quantities {
  array[n] int count_sim;
  for (i in 1:n) {
    int k = categorical_rng(w);
    if (k == 1) count_sim[i] = binomial_rng(total[i], 0.5);
    else if (k == 2) count_sim[i] = beta_binomial_rng(total[i], a, b);
    else count_sim[i] = beta_binomial_rng(total[i], b, a);
  }
}
