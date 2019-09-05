data{
  // training datalist, historical observations 
  int<lower=1> N;         // observations
  int<lower=1> G;         // number of groups
  int<lower=1> K;         // number of individual level predictors
  int<lower=1> J;         // number of group-level predictors
  vector[N] Y;            // all observations
  matrix[N, K] X;         // covariate matrix for fixed effects    
  matrix[N, J] Z;         // covariate matrix for group effects 
  int<lower=1,upper=G> g[N];   // group for individual
}
transformed data{ 
  vector[J] alpha;               // prior on dirichlet distribution
  
  for(i in 1:J)
    alpha[i] = 1; 
}
parameters{
  matrix[J, G] z;
  
  cholesky_factor_corr[J] L_u;    // cholesky factor for correlation for intcpt/slope of year effects
	simplex[J] pi_;                 // pi simplex for diagonal of the covariance matrix 
  real<lower=0> tau;              // scale parameter for covariance matrix


  real<lower=0> sigma;
  vector[K] theta; 
  real<lower=1> nu;          // degrees of freedom

}
transformed parameters{
  matrix[G, J] beta;                   // indiv coeffs by group

  beta = (diag_pre_multiply(pi_*J*tau^2, L_u) * z)';

}
model{
  // Priors
  sigma ~ cauchy(0, 2); 
  theta ~ normal(0, 5);
  nu ~ gamma(2,0.1);
  
  pi_ ~ dirichlet(alpha);       
  tau ~ gamma(1,1);             
  L_u ~ lkj_corr_cholesky(1.0);
  to_vector(z) ~ std_normal(); // unscaled group level effects 

  // Likelihood
  Y ~ student_t_lpdf(nu, X*theta + rows_dot_product(beta[g], Z), sigma);
}
generated quantities{ 
  
  vector[N] Y_hat; 
  
  {
    vector[N] mu_hat; 

    mu_hat = X*theta + rows_dot_product( beta[g], Z); 

    for( i in 1:N)  
      Y_hat[i] = student_t_rng(nu, mu_hat[i], sigma); 
  }
}
