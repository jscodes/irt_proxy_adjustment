//
// This Stan program defines a simple linear regression, with a
// for a single continuous covariate, Z. To be used in PMM. 
//

data {
  int<lower=0> N;
  vector[N] z;
  vector[N] x;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  x ~ normal(alpha + beta * z, sigma);
}



