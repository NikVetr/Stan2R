data {
  int<lower=1> R;
  int<lower=1> C;
  row_vector[C] w;
}
parameters {
  matrix[R, C] Z;
  vector<offset=1, multiplier=2>[R] row_effect;
  real<lower=0> sigma;
  matrix[R, C] Y_raw;  // outcome moved from data, non-centred
}
transformed parameters {
  matrix[R, C] M = rep_matrix(row_effect, C) + Z;
  vector[R] weighted = M * w';
  real grand = sum(M[1]) + sum(M[:, 1]);
  matrix[R, C] Y = M + sigma * Y_raw;
}
model {
  to_vector(Z) ~ std_normal();
  row_effect ~ normal(1, 2);
  sigma ~ exponential(2);
  to_vector(Y_raw) ~ std_normal();
}
