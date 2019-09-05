data{
  // training datalist, historical observations 
  int<lower=0> N;         // observations
  int<lower=0> K;         // number of covariates in the design matrix 
  int<lower=0, upper =1> Z[N]; // 0 (no cover) or 1 (has cover) 
  matrix[N,K] X;          // covariate matrix for fixed effects    
}
parameters{
  // for training data model  
	vector[K] beta;                // fixed effects

}
transformed parameters{
  vector[N] mu;       

  mu = inv_logit(X*beta);
}
model{
  // Priors
  beta ~ normal(0,5);

  // Likelihood
  Z ~ bernoulli_logit(mu);
}
