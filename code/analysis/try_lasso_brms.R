rm(list = ls())
library(future)
library(loo)
library(brms)
library(tidyverse)
source('code/analysis/stan_data_functions.R')

rmse <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  sqrt(mean((yrep_mean - y)^2))
}

testing <- F
n_iter <- 5000 

small <- -0.5
lc <- -0.5
k <- 10
species <- 'PSSP'
vr <- 'growth'
dat <- readRDS('data/temp_data/PSSP_growth_survival_dataframe.RDS')

dat <- 
  dat %>% 
  select( - W ) %>% 
  filter( Period == 'Historical') %>% 
  filter( !is.na(logarea.t1 ), !is.na(logarea.t0)) 

if(testing) { 
  dat <-
    dat %>% 
    group_by(yid) %>% 
    sample_frac(0.25, replace = F) %>% 
    ungroup() 
}

dat$Y <- dat$logarea.t1
dat$size <- dat$logarea.t0
dat$small <- as.numeric(dat$logarea.t0 < small)

dat <- 
  dat %>% 
  mutate( W.total = rowSums( select( . , starts_with('W')))) %>% 
  mutate( W.intra = rowSums( select( . , !!(paste0('W.', species))))) %>% 
  mutate( W.inter = W.total - W.intra)

dat$censored <- 0
dat$censored <- (dat$Y <= lc)*-1 
dat$plant_id <- paste( dat$quad, dat$trackID, sep = '_') 

plot( dat$Y )
abline( h = lc )
plot(dat$size)
abline( h = small )

plot(dat$size, dat$Y)


model_combos <- read_rds('output/model_combos.rds') %>% distinct(model)
C_spec <- paste(model_combos$model[2:8], collapse = '+' )
group_spec <- '(1|yid) + (1|plant_id)'
LHS <- 'Y | cens(censored)'
RHS <- paste( 'size + W.inter + W.intra', C_spec, group_spec, sep = ' + ' )
my_form <- paste(LHS, RHS, sep = ' ~ ')

scaled_covs <- 
  dat %>% 
  select( Y, size, starts_with('C.T'), starts_with('C.V'), starts_with('W.'))

dat <- dat %>% 
  select( yid, quad, trackID, censored ) %>% 
  bind_cols(
    scaled_covs %>% 
    select( Y, size, starts_with('C.T'), starts_with('C.V'), starts_with('W')) %>% 
    mutate_all(scale))

write_rds(dat, path = 'output/brms_kfold/PSSP_growth_data.RDS')

var_names <- names( dat %>% select( size:W.inter) )

my_form <- paste0( LHS, ' ~  (1|yid) + (1|trackID/quad) + ', paste(var_names, collapse = ' + '))

plan(multiprocess)
fit <- brm( my_form, 
     data = dat, 
     family = student, 
     prior = set_prior('horseshoe()', 'b'), 
     future = T, 
     silent = F, 
     control = list(adapt_delta = 0.99, max_treedepth = 20), 
     iter = n_iter)

saveRDS(fit, 'output/brms_kfold/PSSP_horseshoe_fit2.RDS')
