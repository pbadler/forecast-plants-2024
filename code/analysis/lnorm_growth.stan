data{
  // training datalist, historical observations 
  int<lower=0> N;         // observations
  int<lower=0> K;         // number of covariates in the design matrix 
  vector<lower=0>[N] Y;   // all observations
  matrix[N,K] X;          // covariate matrix for fixed effects    

}
parameters{
  // for training data model  
	vector[K] beta;                // fixed effects
	real<lower=0> sigma;
}
transformed parameters{
  vector[N] mu;             

  mu = log(X*beta); 
}
model{
  // Priors
  beta ~ normal(0,5);
  sigma ~ cauchy(0, 2); 
  
  // Likelihood
  Y ~ lognormal(mu, sigma);
}
generated quantities {
  vector[N] Y_hat;
  
  for(i in 1:N) 
    Y_hat[i] = lognormal_rng(mu[i], sigma);
}
