data{
  int<lower=0> N;                   // observations
  int<lower=0> Y[N];                // observation vector
  
  int<lower=0> K;                   // # of fixed effects
  matrix[N,K] X;                    // fixed effects matrix with open space in plot as intercept 
  int<lower=0> G;                   // groups
  int<lower=0> g[N];                // group id

  // Hold data for out of sample prediction 
  int<lower=0> hold_N;                    // observations
  int<lower=0> hold_Y[hold_N];            // observation vector
  matrix[hold_N,K] hold_X;                // fixed effects matrix with open space in plot as intercept 
  int<lower=0> hold_G;                    // groups
  int<lower=0> hold_g[hold_N];            // group id
  
  // IBM data all data for cover predictions 
  int<lower=0, upper = 1> IBM;            // flag whether IBM predictions are generated 
  int<lower=0> IBM_N;                     // observations
  int<lower=0> IBM_Y[IBM_N];              // observation vector
  
  matrix[IBM_N,K] IBM_X;                  // fixed effects matrix with open space in plot as intercept 
  int<lower=0> IBM_G;                     // groups
  int<lower=0> IBM_g[IBM_N];              // group id

}
parameters{
  vector[K] beta;                // fixed effects 
	real<lower=0> tau;             // scale parameter for group-level effects
	vector[G] u_raw;               // raw group-level effects
  real<lower=0> theta;           // negative binomial scale 
}
transformed parameters{
  vector[G] u;                // group-level effects
  vector[N] mu;
  vector[N] fixef;

  for(j in 1:G)
    u[j] = u_raw[j] * tau;
  
  fixef = X*beta;
    
  for(n in 1:N){
    mu[n] = exp(fixef[n] + u[g[n]]);
  }
}
model{
  // Priors
  beta ~ normal(0,1);
  tau ~ cauchy(0,1);                     // gamma as per rstanarm glmer vignette
  u_raw ~ normal(0,1);
  theta ~ cauchy(0,1);

  // Likelihood
  Y ~ neg_binomial_2(mu, theta);
}
generated quantities{
  vector[N] log_lik; 
  vector[N] Y_hat; 
  
  vector[hold_N] hold_mu; 
  vector[hold_N] hold_log_lik; 
  vector[hold_N] hold_Y_hat; 
  vector[hold_N] hold_SE; 
  real hold_SSE; 
  
  vector[IBM_N] IBM_Y_hat; 
  

  for(n in 1:N){ 
    log_lik[n] = neg_binomial_2_lpmf( Y[n] | mu[n], theta); 
    Y_hat[n] = neg_binomial_2_rng( mu[n], theta);
  }
  
  if(hold_N > 2){ 
    vector[hold_N] hold_fixef;   
    vector[hold_G] hold_u; 
  
    for(j in 1:hold_G)
      hold_u[j] = normal_rng(0, 1) * tau;
  
    hold_fixef = hold_X*beta;
    
    for(n in 1:hold_N){
      hold_mu[n] = exp(hold_fixef[n] + hold_u[hold_g[n]]);
      hold_log_lik[n] = neg_binomial_2_lpmf(hold_Y[n] | hold_mu[n], theta); 
      hold_Y_hat[n] = neg_binomial_2_rng(hold_mu[n], theta);
    }
    
  }else if( hold_N <= 2 ){ 
    hold_log_lik = to_vector(rep_array(negative_infinity(), hold_N));
    hold_mu = to_vector(rep_array(0, hold_N));
    hold_Y_hat = to_vector(rep_array(0, hold_N));
  }
  
  for( i in 1:hold_N){
    hold_SE[i] = (hold_mu[i] - hold_Y[i])^2 ;
  }
  hold_SSE = sum(hold_SE)/hold_N;

  if(IBM == 1){ 
    vector[IBM_N] IBM_mu; 
    vector[IBM_N] IBM_fixef;   
    vector[IBM_G] IBM_u; 
  
    for(j in 1:IBM_G)
      IBM_u[j] = normal_rng(0, 1) * tau;
  
    IBM_fixef = IBM_X*beta;
    
    for(n in 1:IBM_N){
      IBM_mu[n] = exp(IBM_fixef[n] + IBM_u[IBM_g[n]]);
      IBM_Y_hat[n] = neg_binomial_2_rng(IBM_mu[n], theta);
    }
  }else if( IBM == 0 ){ 
    IBM_Y_hat = to_vector(rep_array(0, IBM_N));
  }

}

