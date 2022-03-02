rm(list = ls()) 
library(tidyverse)
library(rstan)

rstan_options(auto_write = T)

source('code/analysis/stan_data_functions.R')

datfiles <- dir('data/temp_data/', '*_growth_survival_dataframe.RDS', full.names = T)

my_dat <- readRDS(datfiles[2])

my_dat <- 
  my_dat %>% 
  filter( year < 2000 ) %>%
  filter(!is.na(logarea), !is.na(area.t0))

# Try more complicated model ------------------------------- # 

my_dat$Y <- scale( my_dat$logarea.t1 )
Y <- as.numeric(my_dat$Y)

test <- hist( my_dat$logarea.t1 )
left_cut <- -0.666
abline(v = left_cut, col = 'red', lty = 2)

table( my_dat$logarea.t1[my_dat$logarea.t1 < left_cut  ] )
U <-  max ( Y[ my_dat$logarea.t1 < left_cut ] )

hist( Y )
abline( v = U, col = 'blue', lty = 2)

make_datlist <- function(my_dat, U, formX = ~ size, formZ = ~ size, formE = ~ size ){ 
  cens <- which( Y <= U)
  obs <-  which( Y >  U)

  Y_obs <- Y[obs]
  N_cens <- length(cens)
  N_obs <- length(obs)
  N <- length(Y)
  size <- as.numeric (scale(my_dat$logarea.t0))
  
  my_dat$size <- size 
  
  my_dat$W.intra <-  my_dat$W.ARTR 
  my_dat$W.inter <- rowSums(my_dat[, c('W.HECO', 'W.POSE', 'W.PSSP')])

  X <- model.matrix( formX, data = my_dat)
  Z <- model.matrix( formZ, data = my_dat)
  E <- model.matrix( formE, data = my_dat)

  K <- ncol(X)
  J <- ncol(Z)
  g <- as.numeric( factor( my_dat$year))
  G <- length(unique( g))
  D <- ncol(E)

  pars <- c('N', 'K', 'J', 'N_obs', 'N_cens', 'Y', 'Y_obs', 'X', 'cens', 'obs', 'G', 'g', 'Z', 'U', 'E', 'D')
  datlist <- lapply(pars, FUN = function(x) eval(parse(text = x)))
  names(datlist) <- pars
  return(datlist)
} 


datlist <- make_datlist(my_dat, U, formZ = ~ size)

iter <- 2000
fit1 <- stan('code/analysis/growth2_student_t.stan', 
             data = datlist, 
             cores = 4, 
             iter = iter)

get_elapsed_time(fit1)

mean( colSums( extract( fit1, 'log_lik')$log_lik))  # average log likelihood pointwise probability density 

Y_hat1 <- summary(fit1, 'Y_hat')$summary

plot( Y_hat1[, '50%'], Y)
abline(0,1, col = 'red')

test1 <- test_pp_intervals(Y, fit1)

test1 %>% 
  filter( limit_level == 95) %>% 
  arrange( med ) %>% 
  mutate( id = row_number())  %>% 
  ggplot ( aes( x = id, y = Y)) + 
  geom_errorbar(aes( ymin = lower, ymax = upper), col = 'lightblue')  + 
  geom_point( alpha = 0.1) + 
  geom_line( aes( y = med)) + 
  theme_bw()

test1 %>% 
  filter( Y > U ) %>% 
  group_by( limit_level) %>% 
  summarise( mean ( above ), mean(below), mean( above | below ))

get_year_effects <- function(my_fit, dat, J, G){ 
  u <- summary(my_fit, 'u')$summary[, '50%']
  u <- matrix( u, ncol = J, nrow = G, byrow = T)
   
  u_df <- as.data.frame(u) 
  out <- cbind( data.frame( u_df, year = sort(unique(dat$year))))
  names(out) <- c( paste0('u', 1:J), 'year')
  return(out)
}


ranef1 <- get_year_effects(fit1, my_dat, J = 2, G = datlist1$G)

plot(ranef1$year, ranef1$u1, type = 'l', col = 'red')
points( ranef1$year, ranef1$u2, type = 'l', col = 'blue', lty = 2)
plot(x = ranef1$u1, y = ranef1$u2)

