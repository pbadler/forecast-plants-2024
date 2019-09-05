data{
  // training datalist, historical observations 
  int<lower=0> N;         // observations
  int<lower=0> J;         // number of group level effects
  int<lower=0> K;         // number of covariates in the design matrix 
  vector[N] Y;            // all observations
  row_vector[J] Z[N];     // simple design matrix for random year effects (year x size) 
  int<lower=0> G;         // groups
  int<lower=0> g[N];      // group id
  matrix[N,K] X;          // covariate matrix for fixed effects    
  
  int<lower=0> N_cens;    // censored 
  int<lower=0> N_obs;     // not censored 
  int<lower=0, upper=N> obs[N_obs];  // index of not censored 
  int<lower=0, upper=N> cens[N_cens]; // index of censored observations 
  vector[N_obs] Y_obs;    // only obs with valid size
  real<upper=min(Y_obs) > U;  // Upper limit of censored data 

}
transformed data{ 
  vector[J] alpha;               // prior on dirichlet distribution
  
  for(i in 1:J)
    alpha[i] = 1; 
}
parameters{
  // for training data model  
	vector[K] theta;                // fixed effects
	real<lower=0> sigma;
	real<lower=1> nu;               // degrees of freedom

  cholesky_factor_corr[J] L_u;    // cholesky factor for correlation for intcpt/slope of year effects
	matrix[J,G] u_raw;              // raw group effects 
	simplex[J] pi_;                 // pi simplex for diagonal of the covariance matrix 
  real<lower=0> tau;              // scale parameter for covariance matrix
}
transformed parameters{
  vector[N] ranef;             
  
  {
    vector[J] u[G];   // Group level effects 
    
    // generates covMat diagonal; Explained by rstanarm glmer vignette
    for(j in 1:G)
      u[j] = diag_pre_multiply(pi_*J*tau^2, L_u) * col(u_raw, j);    

    for(n in 1:N)
      ranef[n] = Z[n]*u[g[n]];
  }
  
}
model{
  // Priors
  theta ~ normal(0, 5);
  sigma ~ cauchy(0, 2); 
  nu ~ gamma(2,0.1);

  pi_ ~ dirichlet(alpha);       // dirichlet as per rstanarm glmer vignette
  tau ~ gamma(1,1);             // gamma as per rstanarm glmer vignette
  L_u ~ lkj_corr_cholesky(1.0);
  to_vector(u_raw) ~ std_normal(); // unscaled group level effects 
  
  // Likelihood
  Y_obs ~ student_t_lpdf(nu, (X*theta + ranef)[obs], sigma);
  target += student_t_lcdf( U | nu, (X*theta + ranef)[cens], sigma);  // censored observations 
}
generated quantities{ 
  
  real Y_hat[N] = student_t_rng(nu, X*theta + ranef, sigma);

}
