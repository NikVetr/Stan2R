data {
  int<lower=1> N;
  int<lower=2> K;
  vector<lower=0>[K] alpha;
}
parameters {
  simplex[K] theta;
  simplex[K] phi;
}
model {
  theta ~ dirichlet(alpha);
}
generated quantities {
  array[N] int z_sim;
  for (n in 1:N) z_sim[n] = categorical_rng(theta);
  int z_new = categorical_rng(phi);
}
