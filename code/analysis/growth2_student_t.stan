data{
  // training data
  int<lower=0> N;                     // num. observations
  vector[N] Y;                        // response
  int<lower=1> K;                     // num. of covariates 
  matrix[N,K] X;                      // fixed effects matrix   
  int<lower=1> J;                     // num. of group level effects
  row_vector[J] Z[N];                 // random year effects (year x size)
  int<lower=1> G;                     // num. groups (years)
  int<lower=1> g[N];                  // group id (year id)
  int<lower=0> N_cens;                // censored (Y < U)
  int<lower=0> N_obs;                 // not censored (Y >= U)
  int<lower=0, upper=N> obs[N_obs];   // index of not censored 
  int<lower=0, upper=N> cens[N_cens]; // index of censored observations 
  vector[N_obs] Y_obs;                // obs with valid size
  real<upper=min(Y_obs)> U;           // Upper limit of cens. data 
  
  // holdout data
  int<lower=0> hold_N;              // observations
  vector[hold_N] hold_Y;            // observation vector
  row_vector[J] hold_Z[hold_N];     // simple intercept x size
  int<lower=0> hold_G;              // groups
  int<lower=0> hold_g[hold_N];      // group id
  matrix[hold_N,K] hold_X;          // covariate matrix

  // For generating IBM predictions use all data 
  int<lower=0, upper=1> IBM;        // Flag to generate IBM preds
  int<lower=0> IBM_N;                  
  matrix[IBM_N,K] IBM_X;
  row_vector[J] IBM_Z[IBM_N];
  int<lower=0> IBM_G;
  int<lower=0> IBM_g[IBM_N];
}
transformed data{ 
  vector[J] alpha;               // prior on dirichlet distribution
  
  for(i in 1:J)
    alpha[i] = 1; 
}
parameters{
  // for training data model  
	real<lower=0> sigma;            // prediction error 
	real<lower=1> nu;               // student-t (fat-tails)
  vector[K] beta;                 // fixed effects
  
  cholesky_factor_corr[J] L_u;    // correlation of intcpt/slope 
	matrix[J,G] u_raw;              // scaled group effects 
	simplex[J] pi_;                 // diagonal of the covariance matrix 
  real<lower=0> tau;              // scale for covariance matrix
}
transformed parameters{
  vector[N] mu;             
  vector[J] u[G]; // un-scaled group effects 
  matrix[J,J] Sigma_L = diag_pre_multiply(pi_*J*tau^2, L_u);  // cholesky of covariance 
  
  {
    vector[N] fixef = X*beta;     // fixed effects 

    for(j in 1:G)
      u[j] = Sigma_L * col(u_raw, j);    

    for(n in 1:N)
      mu[n] = fixef[n] + Z[n]*u[g[n]];
  }
  
}
model{
  // Priors
  beta ~ normal(0, 5);
  sigma ~ cauchy(0, 2); 
  nu ~ gamma(2,0.1);
  
  pi_ ~ dirichlet(alpha);       // dirichlet as per rstanarm glmer
  tau ~ gamma(1,1);             // gamma as per rstanarm glmer
  L_u ~ lkj_corr_cholesky(1.0);
  to_vector(u_raw) ~ std_normal();

  // Likelihood
  Y_obs ~ student_t_lpdf(nu, mu[obs], sigma);
  target += student_t_lcdf( U | nu, mu[cens], sigma);  // censored obs.  
}
generated quantities{ 
  real Y_hat[N] = student_t_rng(nu, mu, sigma);
  real log_lik[N]; 
  vector[hold_N] hold_mu;         
  vector[hold_N] hold_log_lik;
  real hold_SSE; 
  vector[IBM_N] IBM_mu;           
  
  for(i in 1:N){ 
    if( Y[i] > U){ 
      log_lik[i] = student_t_lpdf(Y[i] | nu, mu[i], sigma);
    }else if(Y[i] <= U){ 
      log_lik[i] = student_t_lcdf( U | nu, mu[i], sigma);
    }
  }
  
  {
    // Generate out of sample predictions and log-likelihoods 
    vector[J] hold_u[hold_G];       
    vector[hold_N] hold_fixef = hold_X*beta; 
    vector[hold_G*J] ones = rep_vector(1.0, hold_G*J); 
    vector[hold_G*J] zeros = rep_vector(0.0, hold_G*J); 
    matrix[J, hold_G] hold_u_raw = to_matrix( normal_rng(zeros, ones), J, hold_G); 
    vector[hold_N] hold_SE; 

    for(i in 1:hold_G)
      hold_u[i] = Sigma_L * col(hold_u_raw, i);
      
    for(i in 1:hold_N)
      hold_mu[i] = hold_fixef[i] + hold_Z[i]*hold_u[hold_g[i]];
    
    for(i in 1:hold_N){ 
      if( hold_Y[i] > U){ 
        hold_log_lik[i] = student_t_lpdf(hold_Y[i] | nu, hold_mu[i], sigma);
      }else if(hold_Y[i] <= U){ 
        hold_log_lik[i] = student_t_lcdf( U | nu, hold_mu[i], sigma);
      }
      hold_SE[i] = ( hold_mu[i] - hold_Y[i] )^2 ;
    }
    
    hold_SSE = sum(hold_SE);
  }

  if( IBM ==  1 ){ 
    vector[J] IBM_u[IBM_G];       
    vector[IBM_N] IBM_fixef = IBM_X*beta; 
    vector[IBM_G*J] ones = rep_vector(1.0, IBM_G*J); 
    vector[IBM_G*J] zeros = rep_vector(0.0, IBM_G*J); 
    matrix[J, IBM_G] IBM_u_raw = to_matrix( normal_rng(zeros,ones), J, IBM_G); 
  
    for(i in 1:IBM_G)
      IBM_u[i] = Sigma_L * col(IBM_u_raw, i);
      
    for(i in 1:IBM_N)
      IBM_mu[i] = IBM_fixef[i] + IBM_Z[i]*IBM_u[IBM_g[i]];

  }else{
    
    IBM_mu = to_vector(rep_array(0, IBM_N));
  }
}
