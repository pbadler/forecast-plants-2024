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

testing <- T
n_iter <- 2000 

model_combos <- read_rds('output/model_combos.rds')

model_combos <- model_combos %>% distinct(model)
model_combos$C_spec <- model_combos$model
model_combos$group_spec <- '( 1 + size | yid ) + (1 | plant_id)'

group_spec <- c(
  '( 1 + size | yid ) + (1 | plant_id)', 
  '(1 | yid ) + ( 1 | plant_id)', 
  '( 1 + size | yid )', 
  '(1 | yid )')

models <- data.frame( group_spec = group_spec )
models$C_spec <- "NULL"

small <- -0.5
lc <- -0.5
k <- 10
species <- 'ARTR'
vr <- 'growth'
dat <- readRDS('data/temp_data/ARTR_growth_survival_dataframe.RDS')

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
folds <- kfold_split_grouped(k, dat$yid) 
dat$folds <- folds
k_folds <- 
  dat %>% 
  distinct(yid, folds)

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

models <- 
  models %>% 
  mutate( id = row_number())

write_rds(models, path = 'output/brms_kfold/growth_ranef_models.rds')

make_stancode(my_form, 
              data = dat, 
              family  = 'student', 
              prior = set_prior('horseshoe()', 'b') + set_prior('cauchy(0,2)', 'sigma'), 
              iter = n_iter)


for( i in 1:nrow(models)){ 
  
  group_spec <- models$group_spec[i]
  C_spec <- models$C_spec[i]
  LHS <- 'Y | cens(censored)'
  RHS <- 'size + W.inter + W.intra' 
  id <- models$id[i]
  output_file <- paste0('output/brms_kfold/', species, '_', vr, 'ranef_brms_', id, '.RDS')
  
  if( C_spec != 'NULL'){ 
    RHS <- paste( RHS, C_spec, sep = ' + ')
  }
  
  if( group_spec != 'NULL'){ 
    RHS <- paste( RHS, group_spec, sep = ' + ')
  }

  my_form <- as.formula( paste0 (LHS, '~ ', RHS))
  
  
  fit <- brm(my_form, 
             data = dat, 
             family  = 'student', 
             prior = set_prior('normal(0,5)', 'b') + set_prior('cauchy(0,2)', 'sigma'), 
             iter = n_iter)
  
  start <- Sys.time()
  plan(multiprocess)
  kf <- kfold(fit, folds = dat$folds, save_fits = T) 
  et <- Sys.time() - start 
  
  kfp <- kfold_predict(kf)
  error <- rmse(y = kfp$y, yrep = kfp$yrep)
  
  out <- list( id = id, form = my_form, kf = kf, error = error, elapsed_time = et) 
  write_rds( out, path = output_file)
}

