rm(list = ls())
library(pbapply)
library(rstan)
library(tidyverse)
library(loo)

rstan_options(auto_write = T)
options( mc.cores = parallel::detectCores())

source('code/analysis/stan_data_functions.R')

testing <- T

if( testing ){ 
  # STAN pars -------------- 
  n_iter <- 500
  nthin <- 4
  n_mods_per_species <- 2 
}else{
  # STAN pars -------------- 
  n_iter <- 2000
  nthin <- 4
}

# --------------------------
vr <- 'survival'
model_file <- 'code/analysis/survival.stan'
stan_model( file = model_file, model_name = vr)
model_combos <- read_rds('output/model_combos.rds')
load('data/temp_data/prepped_dfs.rda')
total <- length(model_combos)

# Model Parameters 
small <- -1                 ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(paste0 ('~ size + small + W.intra + W.inter + C')) ### Fixed effects design matrix (include climate as "C")

if(testing == T){ 
  # just fit three models per species for testing 
  model_combos <- 
    model_combos %>% 
    filter( species %in% c('ARTR', 'PSSP')) %>% 
    group_by( species) %>% 
    arrange( species, fold, model ) %>% 
    filter( row_number() <= n_mods_per_species ) 
}

model_combos <- 
  model_combos %>% 
  arrange( desc(id)) %>% 
  split(., .$id)

pblapply(model_combos, 
         get_model_score, 
         vr, 
         model_file, 
         prepped_dfs,
         n_iter = n_iter,
         nthin = nthin)



