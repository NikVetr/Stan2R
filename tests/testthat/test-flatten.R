hier <- c("data {", "  int<lower=1> J;", "  vector[J] y;", "}",
          "parameters {", "  real mu;", "  real<lower=0> sd_group, sigma; // two on one line", "  vector[J] z;", "}",
          "model {", "  mu ~ normal(0, 1);", "  sd_group ~ normal(0, 1); // group scale",
          "  for (j in 1:J)", "    z[j] ~ std_normal();", "  sigma ~ exponential(1);",
          "  y ~ normal(mu + sd_group * z, sigma);", "}")

test_that("flatten_model fixes matched scale parameters in transformed data", {
  skip_if_no_stanc()
  out <- strsplit(flatten_model(hier), "\n")[[1]]
  expect_true("transformed data {" %in% out)
  expect_true("  real<lower=0> sd_group = 1;  // set by flatten_model()" %in% out)
  expect_true("  // sd_group ~ normal(0, 1); // group scale" %in% out)
  expect_true("  real<lower=0> sigma; // two on one line" %in% out)
  p <- stan_program(out)
  expect_false("sd_group" %in% p$symbols$name[p$symbols$role == "param"])
})

test_that("flatten_model can move parameters to data and set per-parameter values", {
  skip_if_no_stanc()
  out <- flatten_model(hier, put_params_in_data = TRUE)
  p <- stan_program(out)
  expect_equal(p$symbols$role[p$symbols$name == "sd_group"], "data")
  out2 <- flatten_model(hier, scale_identifier = "^z$|^sd_", set_as = c(sd_group = 2, z = 0))
  expect_match(out2, "vector[J] z = rep_vector(0, J);", fixed = TRUE)
  expect_match(out2, "; // z[j] ~ std_normal();", fixed = TRUE)
  expect_silent(stan_program(out2))
})

test_that("retrieve_block_objects lists declarations and sampled variables", {
  skip_if_no_stanc()
  d <- retrieve_block_objects(hier, "parameters")
  expect_equal(d$var_name, c("mu", "sd_group", "sigma", "z"))
  expect_equal(d$bounds[2], "lower=0, upper=Inf")
  expect_equal(d$dimensions[4], "J")
  s <- retrieve_block_objects(hier, "model", "sampling")
  expect_equal(s$var_name, c("mu", "sd_group", "z", "sigma", "y"))
})
