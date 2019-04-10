data{
  // training datalist, historical observations 
  int<lower=0> N;                       // # observations total 
  int<lower=0, upper=1> S[N];           // success
  int<lower=0> K;                       // # of fixed effects
  matrix[N,K] X;                        // fixed effects matrix   
  int<lower=0> J;                       // # of group level effects
  row_vector[J] Z[N];                   // group level effects matrix
  int<lower=0> G;                       // groups
  int<lower=0> g[N];                    // group id

  // holdout data
  int<lower=0> hold_N;                  // hold_N == 2 if no held out data
  int<lower=0, upper=1> hold_S[hold_N];
  matrix[hold_N,K] hold_X;
  row_vector[J] hold_Z[hold_N];
  int<lower=0> hold_G;
  int<lower=0> hold_g[hold_N];
  
  // For generating IBM predictions use all data 
  int<lower=0, upper=1> IBM;          // Flag if predictions for IBM are to be generated 
  int<lower=0> IBM_N;                  
  matrix[IBM_N,K] IBM_X;
  row_vector[J] IBM_Z[IBM_N];
  int<lower=0> IBM_G;
  int<lower=0> IBM_g[IBM_N];
}
transformed data{ 
  vector[J] a;                          // prior on dirichlet distribution

  for(i in 1:J)
    a[i] = 1; 
}
parameters{
  // for training data model  
	vector[K] beta;                       // fixed effects
	simplex[J] pi_;                       // simplex for diagonal of group-level covariance matrix 
	real<lower=0> tau;                    // scale parameter for group-level covariance matrix
  cholesky_factor_corr[J] L_u;          // cholesky factor for group-level correlation
	matrix[J,G] u_raw;                    // raw group-level effects 
}
transformed parameters{
  // for training data model  
  vector[N] mu;                         // linear predictor 
  vector[J] u[G];                       // group-level effects 
  matrix[J,J] Sigma_L;                  // cholesky of covariance matrix
  vector[J] sigma_j;                    // diagonal of covariance matrix 
    
  sigma_j = pi_*J*tau^2;                      // Explained by rstanarm glmer vignette
  Sigma_L = diag_pre_multiply(sigma_j, L_u);  // multiply variance by correlation  
  
  for(j in 1:G)
    u[j] = Sigma_L * col(u_raw, j);    
  
  {
    vector[N] fixef;              
    fixef = X*beta;
    for(i in 1:N){
      mu[i] = fixef[i] + Z[i]*u[g[i]];
    }    
  }
}
model{
  // Priors
  beta ~ normal(0,5);
  pi_ ~ dirichlet(a);                   // dirichlet as per rstanarm glmer vignette
  tau ~ gamma(1,1);                     // gamma as per rstanarm glmer vignette
  L_u ~ lkj_corr_cholesky(1.0);
  to_vector(u_raw) ~ normal(0,1);

  // Likelihood
  S ~ bernoulli_logit(mu);
}
generated quantities {
  vector[N] log_lik;
  vector[hold_N] hold_log_lik;
  vector[hold_N] hold_mu;           // predicted survival probabilty for held out data
  vector[hold_N] hold_fixef;
  vector[hold_N] hold_SE; 
  real hold_SSE; 
  vector[IBM_N] IBM_mu;     

  for(i in 1:N)
    log_lik[i] = bernoulli_logit_lpmf(S[i] | mu[i]);

  if(hold_N > 2){
    // Run if hold out data is supplied.
    vector[J] hold_u[hold_G];
    matrix[J, hold_G] hold_u_raw;

    for(i in 1:hold_G)
      for(j in 1:J)
        hold_u_raw[j, i] = normal_rng(0,1);

    for(j in 1:hold_G)
      hold_u[j] = Sigma_L * col(hold_u_raw, j);

    hold_fixef = hold_X*beta;

    for(i in 1:hold_N){
      hold_mu[i] = hold_fixef[i] + hold_Z[i]*hold_u[hold_g[i]];
      hold_log_lik[i] = bernoulli_logit_lpmf(hold_S[i] | hold_mu[i]);
    }
    
    }else if(hold_N <= 2 ){
      hold_fixef = to_vector(rep_array(0, hold_N));
      hold_mu = to_vector(rep_array(0, hold_N));
      hold_log_lik = to_vector(rep_array(negative_infinity(), hold_N));
  }
  
  for( i in 1:hold_N){ 
    hold_SE[i] = (inv_logit(hold_mu[i]) - hold_S[i])^2 ;
  }
  
  hold_SSE = sum(hold_SE);
    
  if( IBM ==  1 ){ 

    vector[IBM_N] IBM_fixef;    
    vector[J] IBM_u[IBM_G];
    matrix[J, IBM_G] IBM_u_raw;
  
    for(i in 1:IBM_G)
      for(j in 1:J)
        IBM_u_raw[j, i] = normal_rng(0,1);

    for(j in 1:IBM_G)
      IBM_u[j] = Sigma_L * col(IBM_u_raw, j);

    IBM_fixef = IBM_X*beta;

    for(i in 1:IBM_N){
      IBM_mu[i] = IBM_fixef[i] + IBM_Z[i]*IBM_u[IBM_g[i]];
    }
  }else {
    IBM_mu = to_vector(rep_array(0, IBM_N));
  }
    
}
