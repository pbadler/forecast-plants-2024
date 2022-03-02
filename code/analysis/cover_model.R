rm(list = ls())
library(tidyverse)
library(rstan)
cover <- readRDS('data/temp_data/all_cover.RDS')
clim <- readRDS('data/temp_data/all_clim_covs.RDS')
quads <- read.csv('data/quad_info.csv')

logit <- function( p ) { log ( p/(1-p)) }
inv_logit <- function(x) {1/(1 + exp(-x))}

make_cov_data_list <- function(cover, climate, formX, formZ){ 
  
  x <- 
    cover %>% 
      left_join(cover %>% mutate(year = year + 1, cover0 = cover ) %>% 
                select(-cover, -area) ) %>% 
      arrange( QuadName, spp) %>% 
      select( -area ) %>% 
      left_join(clim, by = c('year', 'Treatment')) %>% 
      filter( !is.na(cover), !is.na(cover0) )
    
  X <- model.matrix(data = x, formX )
  Y_obs <- x$cover[x$cover > 0 ]
  Y <- x$cover  
  P <- as.numeric( x$cover > 0 )
  obs <- which(x$cover > 0)
  N <- length(P)
  N_obs <- length(Y_obs)
  K <- ncol(X)
  
  year <- x$year
  Treatment <- x$Treatment 
  q <- x$QuadName 
  
  Z <- model.matrix(data = x, formZ)
  g <- as.numeric( factor(x$year))
  J <- ncol(Z)
  G <- length(unique(g))
  
  pars <- c('X', 'Y', 'P', 'Y_obs', 'N_obs', 'obs', 'N', 'K', 'Z', 'g', 'J', 'G', 'q', 'year', 'Treatment')
  
  out <- lapply( pars, function(x) eval(parse(text = x))) 
  names(out ) <- pars 
  
  return(out )
}

get_year_effects <- function(my_fit, dat, re_name = 'u', J = 1){ 
  us <- summary(my_fit, 'u')$summary
  pats <- paste0( 1:J, ']')
  u <- matrix(NA, ncol = J, nrow = nrow(us))
  for( i in 1:J){ 
    u[,i] <- us [ row.names(us) [ grep ( pattern = pats[i], row.names(us) ) ] , '50%']
  }
  u_df <- as.data.frame(u) 
  out <- cbind( data.frame( u_df, year = sort(unique(dat$year))))
  names(out) <- c( paste0('u', 1:J), 'year')
  return(out)
}
cover <- 
  cover %>% 
  filter( !is.na(cover) ) %>%
  split(.$spp)

dat4stan <- make_cov_data_list(cover[[2]], climate, cover ~ cover0, cover ~ cover0 )
dat4stan2 <- make_cov_data_list(cover[[2]], climate, cover ~ cover0, cover ~ 1)

fit2 <- stan('code/analysis/cover_model.stan', data = dat4stan, cores = 4, chains = 4)

fit3 <- stan('code/analysis/cover_model.stan', data = dat4stan2, cores = 4, chains = 4)

summary(fit2, 'beta')$summary
Y_hat <- summary(fit2, 'Y_hat')$summary[, '50%']
Y_hat2 <- summary(fit3, 'Y_hat')$summary[, '50%']

sqrt( mean( (Y_hat - dat4stan$Y)^2))
sqrt( mean( (Y_hat2 - dat4stan$Y)^2))

plot(Y_hat, dat4stan$Y)
abline(0, 1)

my_fit <- fit2
dat <- dat4stan2
J <- 2



end_year <- 2017
RE_est1 <- get_year_effects(fit2, dat = dat4stan, re_name = 'u', J = 2)
RE_est2 <- get_year_effects(fit3, dat = dat4stan2, J = 1)

data.frame(plot = dat4stan$q, cover = dat4stan$Y, cover0 = as.numeric (dat4stan$X[, 2]), year = dat4stan$year) %>% 
  mutate( growth = log( cover ) - log( cover0 ) ) %>% 
  filter( is.finite(growth) ) %>% 
  filter( year < end_year ) %>% 
  ggplot( aes( x = year, y = growth )) + 
  geom_point() + 
  stat_summary(fun.y = 'mean', geom = 'line', color = 'red') + 
  geom_line( data = RE_est1 %>% filter( year < end_year), 
             aes( x = year, y = u1), col = 'blue') + 
  geom_line( data = RE_est1 %>% filter( year < end_year), 
             aes( x = year, y = 10*u2), col = 'pink') + 
  geom_line( data = RE_est2 %>% filter( year < end_year), 
             aes( x = year, y= u1), col = 'green')


ggplot(data = RE_est1 %>% filter( year < end_year), aes( x = year, y = u2)) + geom_line()

dat4stanzeros <- make_cov_data_list(my_cov_zeros, cover ~ cover0)

fitZeros <- stan('code/analysis/has_plants.stan', data = dat4stanzeros$ARTR, chains = 4)
summary(fitZeros, 'beta')$summary
mu_hat<-summary(fitZeros, 'mu')$summary 
mu_hat <- mu_hat[, '50%']

betas <- summary(fitZeros, 'beta')$summary[, '50%']

mu <- inv_logit( dat4stan$ARTR$X %*% betas )

length(mu)
length(Y_hat)

plot(Y_hat, dat4stan$ARTR$Y)
points( mu*Y_hat, dat4stan$ARTR$Y, col = 'red', pch = 1)
abline(0,1)

plot(Y_hat[Y_hat < 1 ], dat4stan$ARTR$Y[ Y_hat < 1 ])
points( (mu*Y_hat)[ Y_hat < 1 ], dat4stan$ARTR$Y[Y_hat < 1], col = 'red', pch = 1)
abline(0,1)

plot( mu_hat, dat4stanzeros$ARTR$Z)


curve(logit, 0, 1)
curve(inv_logit, -10, 10)

