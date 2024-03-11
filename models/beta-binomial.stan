data {
  int<lower=0> n;
  array[n] int<lower=0> k;
  array[n] int<lower=0, upper=k> x;
}

transformed data {
  vector[n] log_x = to_vector(log(x));
}

parameters {
  real theta_intercept;
  real theta_slope;
  real<lower=0> conc_offset;
}

transformed parameters {
  vector[n] theta = theta_intercept + rep_vector(theta_slope, n) .* log_x;
  vector<lower=0, upper=1>[n] mu = inv_logit(theta);
  real<lower=2> conc = conc_offset + 2;
  vector<lower=0>[n] shape1 = mu .* rep_vector(conc, n);
  vector<lower=0>[n] shape2 = (1-mu) .* rep_vector(conc, n);
}

model {
  conc_offset ~ exponential(0.1);
  x ~ binomial(k, shape1, shape2);
}
