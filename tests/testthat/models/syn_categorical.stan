// simplex prior and categorical outcomes; the simplex has no explicit prior in one place
data {
  int<lower=1> N;
  int<lower=2> K;
  array[N] int<lower=1, upper=K> z;
  vector<lower=0>[K] alpha;
}
parameters {
  simplex[K] theta;
  simplex[K] phi;
}
model {
  theta ~ dirichlet(alpha);
  z ~ categorical(theta);
}
generated quantities {
  int z_new = categorical_rng(phi);
}
