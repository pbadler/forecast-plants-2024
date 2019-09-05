data{
  // training datalist, historical observations 
  int<lower=0> N;         // All observations
  int<lower=0> N_obs;     // Num. with cover > 0  
  int<lower=0, upper=N> obs[N_obs];   // index of those cover > 0
  vector[N_obs] Y_obs;        // Obs with cover > 0
  int<lower=0, upper=1> P[N]; // 0,1 indicated whether cover is greater than zero 
  int<lower=0> K;         // number of covariates in the design matrix 
  matrix[N,K] X;          // covariate matrix for fixed effects    
  
  // Group level effects   
  int<lower=0> J;         // number of group level effects
  row_vector[J] Z[N];     // simple design matrix for random year effects (year x cover) 
  int<lower=0> G;         // groups
  int<lower=0> g[N];      // group id
  
  
  // Validation data, held out 
  // int<lower=0> hold_N;         
  // int<lower=0> hold_N_obs;     
  // int<lower=0, upper=hold_N> hold_obs[hold_N_obs];   
  // vector[hold_N_obs] hold_Y_obs;        
  // int<lower=0, upper=1> hold_P[hold_N]; 
  // matrix[hold_N,K] hold_X;          
  
}
transformed data{ 
  vector[J] alpha;               // prior on dirichlet distribution
  
  for(i in 1:J)
    alpha[i] = 1; 
}
parameters{
  // for training data model  
	vector[K] beta;                 // fixed effects
	real<lower=0> sigma;            // dispersion parameter
  vector[K] theta;                // fixed effects for binomial model 
  
  cholesky_factor_corr[J] L_u;    // cholesky factor for correlation for intcpt/slope of year effects
	matrix[J,G] u_raw;              // raw group effects 
	simplex[J] pi_;                 // pi simplex for diagonal of the covariance matrix 
	real<lower=0> tau;              // scale parameter for covariance matrix

}
transformed parameters{
  vector[N] mu;                   //  all observations 
  vector[N_obs] mu_obs;           // only cover > 0 
  vector[N] phi;                  //  cover 
  
  vector[J] u[G];                 // scaled and correlated group effects 
  matrix[J,J] Sigma_L;            // cholesky of covariance matrix
  
  vector[N] fixef;                  
  
  fixef = log(X*beta);
  phi = X*theta; 
  
  Sigma_L = diag_pre_multiply(pi_*J*tau^2, L_u);   // generates covMat diagonal; Explained by rstanarm glmer vignette

  for(j in 1:G)
    u[j] = Sigma_L * col(u_raw, j);    
      
  
  for(i in 1:N_obs)   
    mu_obs[i] = fixef[obs[i]] + Z[obs[i]]*u[g[obs[i]]];
  
  for(i in 1:N)   
    mu[i] = fixef[i] + Z[i]*u[g[i]]; 
  
  
}
model{
  // Priors
  beta ~ normal(0,5);
  sigma ~ cauchy(0,2); 
  theta ~ normal(0,5);
  pi_ ~ dirichlet(alpha);       // dirichlet as per rstanarm glmer vignette
  tau ~ gamma(1,1);             // gamma as per rstanarm glmer vignette
  L_u ~ lkj_corr_cholesky(1.0);
  to_vector(u_raw) ~ normal(0,1);

  // Likelihood
  Y_obs ~ lognormal( mu_obs, sigma ); 
  P ~ bernoulli_logit( phi ); 
}
generated quantities{ 
  vector[N] Y_hat;

  for( i in 1:N )
    Y_hat[i] = inv_logit(phi[i])*lognormal_rng(mu[i], sigma);    // Cover taking into account potential loss of all plants 
  
}
