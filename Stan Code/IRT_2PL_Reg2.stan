
// 2-Parameter Logistic IRT model fit on all patient data
// with a single continuous covariate and a categorical 
// indicator for response mode mechanism

data {
  int<lower=1> Obs;               // number of patients by item observations
  int<lower=1> N;                 // number of patients
  int<lower=1> P;                 // number of items
  int<lower=1, upper=N> ii[Obs];  // patient i for obs
  int<lower=1, upper=P> kk[Obs];  // item k for obs
  vector[N] z;                    // subject-level covariate
  vector[N] w;                    // subject-level response mode mechanism
  int<lower=0, upper=1> x[Obs];   // outcome value for patient observations

}
parameters {
  vector[N] theta;
  vector[P] beta;
  vector<lower=0>[P] alpha;     
  real<lower=0> sigma_beta;
  real<lower=0> sigma_alpha;
  real delta;
  real rho;
}
model {
  theta ~ normal(delta*z + rho*w, 1);
  beta ~ normal(0, sigma_beta);
  alpha ~lognormal(0, sigma_alpha);
  sigma_beta ~ cauchy(0, 5);
  sigma_alpha ~ cauchy(0, 5);  
  x ~ bernoulli_logit(alpha[kk] .* (theta[ii] - beta[kk]));
}



