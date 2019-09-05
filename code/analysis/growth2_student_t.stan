data{
  // training data
  int<lower=0> N;                     // observations
  vector[N] Y;                        // all observations
  int<lower=1> K;                     // number of covariates in the design matrix 
  matrix[N,K] X;                      // fixed effects matrix   
  int<lower=1> J;                     // number of group level effects
  row_vector[J] Z[N];                 // simple design matrix for random year effects (year x size) 
  int<lower=1> G;                     // num. groups (years)
  int<lower=1> g[N];                  // group id (year id)
  int<lower=0> N_cens;                // censored 
  int<lower=0> N_obs;                 // not censored 
  int<lower=0, upper=N> obs[N_obs];  // index of not censored 
  int<lower=0, upper=N> cens[N_cens]; // index of censored observations 
  vector[N_obs] Y_obs;    // only obs with valid size
  real<upper=min(Y_obs)> U;  // Upper limit of censored data 
  
  
}
transformed data{ 
  vector[J] alpha;               // prior on dirichlet distribution
  
  for(i in 1:J)
    alpha[i] = 1; 
}
parameters{
  // for training data model  
	real<lower=0> sigma;            // prediction error 
	real<lower=1> nu;               // student-t degrees of freedom (fat-tails)
  vector[K] beta;                 // fixed effects
  
  cholesky_factor_corr[J] L_u;    // correlation of intcpt/slope of group effects
	matrix[J,G] u_raw;              // scaled group effects 
	simplex[J] pi_;                 // diagonal of the covariance matrix 
  real<lower=0> tau;              // scale parameter for covariance matrix
}
transformed parameters{
  vector[N] mu;             
  vector[J] u[G];                 // un-scaled group effects 
 
  {
    vector[N] fixef = X*beta;     // fixed effects 

    // generates covMat diagonal; Explained by rstanarm glmer vignette
    matrix[J,J] Sigma_L = diag_pre_multiply(pi_*J*tau^2, L_u);  // cholesky of covariance matrix
  
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
  
  pi_ ~ dirichlet(alpha);       // dirichlet as per rstanarm glmer vignette
  tau ~ gamma(1,1);             // gamma as per rstanarm glmer vignette
  L_u ~ lkj_corr_cholesky(1.0);
  to_vector(u_raw) ~ std_normal();

  // Likelihood
  Y_obs ~ student_t_lpdf(nu, mu[obs], sigma);
  target += student_t_lcdf( U | nu, mu[cens], sigma);  // integrate out the censored observations 

}
generated quantities{ 
  
  real Y_hat[N] = student_t_rng(nu, mu, sigma);
  real log_lik[N]; 

  for(i in 1:N){ 
    
    if( Y[i] > U){ 
      log_lik[i] = student_t_lpdf(Y[i] | nu, mu[i], sigma);
      
    }else if(Y[i] <= U){ 
      log_lik[i] = student_t_lcdf( U | nu, mu[i], sigma);
    }
  }
  
}
