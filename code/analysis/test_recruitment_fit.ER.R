library(tidyverse)
library(rstan)

load('data/temp_data/prepped_dfs.rda')

dat <- 
  prepped_rec_dfs$HECO %>% 
  filter( year < 2000 )

source('code/analysis/stan_data_functions.R')

mc <- readRDS('output/model_combos.rds')

formC <- as.formula(paste0 ( '~ - 1 + ', mc$model[1]))
formX <- as.formula(' ~ C ')
formZ <- as.formula(' ~ 1')


dl <- process_data(dat, formZ, formX, formC, IBM = 0, vr = 'recruitment', hold = 0)

library(rstanarm)

my_dat <- as.data.frame( dl[c('Y', 'X', 'Z', 'g')] )
my_dat$g


nb_fit <- stan_glmer.nb('Y ~ (1|g)', data = my_dat, cores = 4)

test <- nb_fit$stanfit

test2 <- extract(test)

hist( test2$aux )

extract(test, 'b')

?stan_glm.nb
stan_diag( nb_fit )
stan_trace(nb_fit)

Y_hat <- posterior_predict(nb_fit)
med <- apply( Y_hat, 2, quantile, 0.5 )
low <- apply( Y_hat, 2, quantile, 0.25)
hi <- apply( Y_hat, 2, quantile, 0.75)
Y_hat
N <- ncol(Y_hat)
x <- 1:N

plot( x, med, ylim = c(0, 30), type = 'n')
segments(x0 = x, y0 = low, x1 = x, y1 = hi, col = 'blue', lty = 1)
points( x, dl$Y)


nb_fit$coefficients
nb_fit$stan_summary
nb_fit$family
nb_fit$coefficients

