data {
  int<lower=0> n;
  array[n] int<lower=1> k;
  array[n] int<lower=0, upper=k> x;
  int<lower=1> n_groups;
  array[n] int<lower=1, upper=n_groups> group;
}

transformed data {
  vector<lower=0>[n] k_scaled = to_vector(k) ./ (max(k) * 1.0);
}

parameters {
  real theta_intercept;
  real theta_slope;
  vector[n_groups] theta_group; //this is a test comment
  real<lower=0> theta_group_sd;
  real<lower=0> conc_offset;
}

transformed parameters {
  vector[n] theta = theta_intercept + rep_vector(theta_slope, n) .* k_scaled + theta_group[group] * theta_group_sd;
  vector<lower=0, upper=1>[n] mu = inv_logit(theta);
  real<lower=2> conc = conc_offset * 3 + 40;
  vector<lower=0>[n] shape1 = 1 + mu .* rep_vector(conc, n); //this is another test comment
  vector<lower=0>[n] shape2 = 1 + (1-mu) .* rep_vector(conc, n);
}

model {
  //priors
  theta_intercept ~ std_normal();
  theta_slope ~ std_normal();
  theta_group ~ std_normal();
  theta_group_sd ~ std_normal();
  conc_offset ~ exponential(0.1);
  
  //likelihood
  x ~ beta_binomial(k, shape1, shape2);
}
