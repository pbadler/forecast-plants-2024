data{
  // training data
  int<lower=0> N;                     // observations
  vector[N] Y;                        // all observations
  int<lower=0, upper = 1> P[N];       // presence/absence (1, 0)
  int<lower=1> K;                     // number of covariates
  int<lower=1> L;                     // presence absence covariates. 
  matrix[N,K] X;                      // covariates matrix for cover
  matrix[N,L] W;                      // covariates matrix for pres./abs.
  int<lower=1> J;                     // number of group level effects
  row_vector[J] Z[N];                 // simple design matrix for random year effects (year x size) 
  int<lower=1> G;                     // num. groups (years)
  int<lower=1> g[N];                  // group id (year id)
  int<lower=0> N_cens;                // censored 
  int<lower=0> N_obs;                 // positive abundance 
  int<lower=0, upper=N> obs[N_obs];  // index of not censored 
  int<lower=0, upper=N> cens[N_cens]; // index of censored observations 
  vector<lower=0> [N_obs] Y_obs;    // only obs with valid size

  // holdout data
  int<lower=0> hold_N;              // observations
  int<lower=0> hold_N_obs;          // positive abundance 
  vector[hold_N] hold_Y;            // observation vector
  int<lower=0, upper = 1> hold_P[N];// presence/absence (1, 0)
  row_vector[J] hold_Z[hold_N];     // simple intercept x size design   
  int<lower=0> hold_G;              // groups
  int<lower=0> hold_g[hold_N];      // group id
  matrix[hold_N,K] hold_X;          // covariate matrix
  matrix[hold_N,L] hold_W;          // covariates matrix for pres./abs.
  int<lower=0, upper=hold_N> hold_obs[hold_N_obs];  // index of not censored 
  vector<lower=0> [hold_N_obs] hold_Y_obs; // only obs with valid size

}
transformed data{ 
  vector[J] alpha;               // prior on dirichlet distribution
  
  for(i in 1:J)
    alpha[i] = 1; 
}
parameters{
	vector[K] beta;                 // fixed effects
	real<lower=0> sigma;            // dispersion parameter
  vector[L] theta;                // fixed effects for binomial model 
  
  cholesky_factor_corr[J] L_u;    // correlation for intcpt/slope
	matrix[J,G] u_raw;              // raw group effects 
	simplex[J] pi_;                 // scale correlation 
	real<lower=0> tau;              // scale parameter for covariance matrix
}
transformed parameters{
  vector[N] mu;             
  vector[N] phi = W*theta; 
  vector[J] u[G]; 
  matrix[J,J] Sigma_L = diag_pre_multiply(pi_*J*tau^2, L_u);  
  
  {
    vector[N] fixef = log(X*beta); 
    
    for(j in 1:G)
      u[j] = Sigma_L * col(u_raw, j);    

    for(n in 1:N)
      mu[n] = fixef[n] + Z[n]*u[g[n]];
  }
}
model{
  // Priors
  beta ~ normal(0,5);
  sigma ~ cauchy(0,2); 
  theta ~ normal(0,5);
  pi_ ~ dirichlet(alpha);       
  tau ~ gamma(1,1);             
  L_u ~ lkj_corr_cholesky(1.0);
  to_vector(u_raw) ~ std_normal();

  // Likelihood
  Y_obs ~ lognormal( mu[obs], sigma ); 
  P ~ bernoulli_logit_lpmf( phi ); 
}
generated quantities{ 
  vector[N] Y_hat;
  vector[N_obs] log_lik1; 
  vector[N] log_lik2; 
  
  // Generate out of sample predictions and log-likelihoods 
  vector[hold_N] hold_mu;
  vector[hold_N] hold_phi = hold_W*theta; 
  vector[hold_N_obs] hold_log_lik1; 
  vector[hold_N] hold_log_lik2;
  real hold_SSE; 

  for( i in 1:N ){ 
    Y_hat[i] = inv_logit(phi[i])*lognormal_rng(mu[i], sigma); // preds. 
    log_lik2[i] = bernoulli_logit_lpmf( P[i] | phi[i]);
  }    
  
  for( i in 1:N_obs)
    log_lik1[i] = lognormal_lpdf( Y_obs[i] | mu[obs[i]], sigma ); 
  
  {
    vector[J] hold_u[hold_G];       
    vector[hold_N] hold_fixef = log( hold_X*beta ); 
    vector[hold_G*J] ones = rep_vector(1.0, hold_G*J); 
    vector[hold_G*J] zeros = rep_vector(0.0, hold_G*J); 
    matrix[J, hold_G] hold_u_raw = to_matrix( normal_rng(zeros, ones), J, hold_G); 
    vector[hold_N] hold_SE; 

    for(i in 1:hold_G)
      hold_u[i] = Sigma_L * col(hold_u_raw, i);
      
    for(i in 1:hold_N){ 
      hold_mu[i] = hold_fixef[i] + hold_Z[i]*hold_u[hold_g[i]];
      hold_log_lik2[i] = bernoulli_logit_lpmf(hold_P[i] | hold_phi[i]);
      hold_SE[i] = (inv_logit(hold_phi[i])*exp(hold_mu[i]) - hold_Y[i] )^2;
    }
    
    for(i in 1:hold_N_obs)
      hold_log_lik1[i] = lognormal_lpdf(hold_Y_obs[i] | hold_mu[hold_obs[i]], sigma);
      
    hold_SSE = sum(hold_SE);
  }

  
}
