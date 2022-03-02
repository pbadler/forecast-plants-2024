rm(list = ls())
library(rstan)

schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))


fit <- stan(file = 'code/test/8schools.stan', data = schools_dat, refresh = 0, cores = 4, iter = 500)
id <- 3
nd <- get_num_divergent( fit)
if( nd > 0 ) stop( paste0( nd,  ' Divergent transitions in model run ', id ))

