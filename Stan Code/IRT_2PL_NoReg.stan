
// Modified from Fitting Bayesian item response models in Stata and Stan, 
// Grant et al.

// 2-Parameter IRT model with no covariates
data {
  int<lower=1> Obs;               // number of patients by item observations
  int<lower=1> N;                 // number of patients
  int<lower=1> P;                 // number of items
  int<lower=1, upper=N> ii[Obs];  // patient i
  int<lower=1, upper=P> kk[Obs];  // item k
  int<lower=0, upper=1> x[Obs];   // outcome value for patient observations
}
parameters {
  vector[N] theta;
  vector[P] beta;
  vector<lower=0>[P] alpha;     
  real<lower=0> sigma_beta;
  real<lower=0> sigma_alpha;
}
model {
  theta ~ normal(0, 1);
  beta ~ normal(0, sigma_beta);
  alpha ~lognormal(0, sigma_alpha);
  sigma_beta ~ cauchy(0, 5);
  sigma_alpha ~ cauchy(0, 5);  
  x ~ bernoulli_logit(alpha[kk] .* (theta[ii] - beta[kk]));
}


