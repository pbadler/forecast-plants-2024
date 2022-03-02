functions{
  real age_curve(real x, real a, real b, real c, real min_size){
    return a*exp(-(x-b)^2/(2*c^2)) + min_size;
  }
}
data{
  // training datalist, historical observations 
  int<lower=0> N;         // observations
  int<lower=0> P; 
  // int<lower=0> K;         // number of covariates in the design matrix 
  vector[N] Y;            // all observations
  vector<lower=1>[N] age; 
  int<lower=1, upper=P> p[N]; 
  real<lower=0> min_size; 
  
  // int<lower=0> J;         // number of group level effects
  // row_vector[J] Z[N];     // simple design matrix for random year effects (year x size) 
  // int<lower=0> G;         // groups
  // int<lower=0> g[N];      // group id
}
transformed data{ 
  // vector[J] alpha;               // prior on dirichlet distribution
  // 
  // for(i in 1:J)
  //   alpha[i] = 1; 
}
parameters{
  real<lower=0> a;
  real b; 
  real<lower=0> c; 
  real<lower=0> sigma;
//   cholesky_factor_corr[J] L_u;    // cholesky factor for correlation for intcpt/slope of year effects
// 	matrix[J,G] u_raw;              // raw group effects 
// 	simplex[J] pi_;                 // pi simplex for diagonal of the covariance matrix 
// 	real<lower=0> tau;              // scale parameter for covariance matrix
}
transformed parameters{
  real mu[N];             
  // vector[J] u[G];                 // scaled and correlated group effects 
  // matrix[J,J] Sigma_L;            // cholesky of covariance matrix
  // vector[J] sigma_j;   
  // 
  // sigma_j = pi_*J*tau^2;        // generates covMat diagonal; Explained by rstanarm glmer vignette
  // Sigma_L = diag_pre_multiply(sigma_j, L_u);  // multiply variance by correlation  
  // 
  // for(j in 1:G)
  //   u[j] = Sigma_L * col(u_raw, j);    

  for( n in 1:N)
    mu[n] = log(age_curve(age[n], a, b, c, min_size)); 
  
}
model{

  // Priors
  a ~ cauchy(0,2); 
  b ~ normal(0,5);
  c ~ gamma(1, 1); 
  sigma ~ cauchy(0,2); 
  
  // pi_ ~ dirichlet(alpha);       // dirichlet as per rstanarm glmer vignette
  // tau ~ gamma(1,1);             // gamma as per rstanarm glmer vignette
  // L_u ~ lkj_corr_cholesky(1.0);
  // to_vector(u_raw) ~ normal(0,1);
		
  // Likelihood
  Y ~ normal(mu, sigma);
}
