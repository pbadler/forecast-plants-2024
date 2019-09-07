rm(list = ls())

library(rstan)
library(tidyverse)
library(loo)

rstan_options(auto_write = T)
options( mc.cores = parallel::detectCores())

source('code/analysis/stan_data_functions.R')

testing <- F

if( testing ){ 
  # STAN pars -------------- 
  n_iter <- 500
  nthin <- 4
  
}else{
  # STAN pars -------------- 
  n_iter <- 2000
  nthin <- 4
}

# --------------------------
vr <- 'growth'
model_file <- 'code/analysis/growth2_student_t.stan'
stan_model( file = model_file, model_name = vr, save_dso = T)
model_combos <- read_rds('output/model_combos.rds')
load('data/temp_data/prepped_dfs.rda')
total <- length(model_combos)

# Model Parameters 
small <- -1                 ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(paste0 ('~ size + small + W.intra + W.inter + C')) ### Fixed effects design matrix (include climate as "C")

model_combos <- model_combos %>% split(., .$id)

test <- model_combos[185:length(model_combos)]

pblapply(test, 
         get_model_score, 
         vr, 
         model_file, 
         prepped_dfs,
         n_iter = n_iter,
         nthin = nthin)

