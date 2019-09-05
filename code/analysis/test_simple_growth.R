rm(list = ls()) 
library(tidyverse)
library(rstan)

rstan_options(auto_write = T, cores = 4)
source('code/analysis/stan_data_functions.R')

datfiles <- dir('data/temp_data/', '*_growth_survival_dataframe.RDS', full.names = T)

make_testing_data <- function(datfilename, left_censor, formX = ~ size, formZ = ~ size ){  
  
  my_dat <- readRDS(datfilename)

  my_dat <- 
    my_dat %>% 
    filter( year < 2000, logarea.t1 > left_censor ) %>%
    filter(!is.na(logarea.t1), !is.na(logarea.t0))
  
  my_dat$Y <- scale( my_dat$logarea.t1 )
  Y <- as.numeric(my_dat$Y)

  N <- length(Y)
  my_dat$size <- as.numeric (scale(my_dat$logarea.t0))
  
  X <- model.matrix( formX, data = my_dat)
  Z <- model.matrix( formZ, data = my_dat)

  K <- ncol(X)
  J <- ncol(Z)
  g <- as.numeric( factor( my_dat$year))
  G <- length(unique( g))

  pars <- c('N', 'K', 'J', 'Y', 'X', 'G', 'g', 'Z' )
  datlist <- lapply(pars, FUN = function(x) eval(parse(text = x)))
  names(datlist) <- pars
  return(datlist)
} 

my_dat <- readRDS(datfiles[3])

#my_dat %>% group_by(logarea.t1) %>% summarise( n=  n() )  %>% arrange( desc( n) ) %>% head( 20 )

datlist <- make_testing_data(datfilename = datfiles[2], left_censor = -1, formZ = ~ size)

iter <- 2000

fit1 <- stan('code/analysis/growth2_student_t_simple.stan', data = datlist, cores = 4, iter = iter)
fit2 <- stan('code/analysis/growth2_student_t_optim.stan', data = datlist, cores =4, iter = iter)

get_elapsed_time(fit1)
get_elapsed_time(fit2)

beta <- summary( fit2, 'beta')$summary[, 1]
beta <- matrix( beta, length( beta )/2, 2, byrow = T)
theta <- summary(fit2, 'theta')$summary[,1]
Y_hat2 <- summary(fit2, 'Y_hat')$summary[,1]
Y_hat1 <- summary(fit1, 'Y_hat')$summary[,1]
plot( datlist$Y, Y_hat2)
points(datlist$Y, Y_hat1, col = 'blue')
abline(0,1, col = 'red')

pi_ <- summary(fit, 'pi_')$summary[, 1] 
L_u <- summary(fit, 'L_u')$summary[, 1]
tau <- summary(fit, 'tau')$summary[, 1]
L_u <- matrix( L_u, 2, 2, byrow = T)
u_raw <- summary(fit, 'u_raw')$summary[, 1]

cbind( u_raw[ which( seq_along(u_raw) %% 2 == 1  )] , u_raw[ which( seq_along(u_raw) %% 2 == 0 )])

u_raw <- matrix( u_raw, nrow = 2, ncol = datlist$G , byrow = T) 

t(u_raw)
alg_dat <- list( J = datlist$J, G = datlist$G, pi_ = pi_, L_u = L_u, tau = tau, u_raw = t(u_raw))

test <- stan('code/analysis/stan_lin_alg_operation.stan', data = alg_dat,algorithm = 'Fixed_param')

summary(test, 'Sigma_L')$summary[, 1]
summary(fit, 'Sigma_L')$summary [, 1]

u1 <- matrix( summary(fit, 'u')$summary[, 1], datlist$J, datlist$G, byrow = T)
u2 <- matrix( summary(test, 'u')$summary[, 1], datlist$J, datlist$G, byrow = T)

cbind( t(u1), t(u2))
plot( u1, u2)
