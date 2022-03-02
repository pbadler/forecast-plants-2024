rm(list = ls())
library(pbapply)
library(rstan)
library(tidyverse)
library(loo)
library(rstanarm)

rstan_options(auto_write = T)
options( mc.cores = parallel::detectCores())


source('code/analysis/stan_data_functions.R')

testing <- T
k <- 10 

if( testing ){ 
  # STAN pars -------------- 
  n_iter <- 1000
  nthin <- 4
  n_mods_per_species <- 2 
}else{
  # STAN pars -------------- 
  n_iter <- 2000
  nthin <- 4
}

# --------------------------
vr <- 'recruitment'
#model_file <- 'code/analysis/recruitment.stan'
#stan_model( file = model_file, model_name = vr)

dat <- readRDS('data/temp_data/ARTR_recruitment_dataframe.RDS')

train <- dat[ dat$Period == 'Historical',  ] 

folds <- kfold_split_grouped(k, train$yid) 
train$folds <- folds

k_folds <- 
  train %>% 
  distinct(yid, folds)

train <- 
  train %>% 
  ungroup %>% 
  rowwise( ) %>% 
  mutate( total_basal_cover = cov.HECO + cov.POSE + cov.PSSP) %>% 
  mutate( total_open = 100*100 - total_basal_cover) %>% 
  mutate( l_open = log(total_open))

train$W <- scale( train$total_basal_cover)

train <- as.data.frame( train )   

fit <- stan_glmer.nb( Y ~ W + (1|yid), data = train ) 
fit2 <- stan_glmer.nb(Y ~ W + C.T.1*C.VWC.1 + (1|yid), data = train)

k_res1 <- kfold(fit, K = 10, folds = folds)
k_res2 <- kfold(fit2, K = 10, folds = folds)

compare_models(k_res1, k_res2, detail = T)


data.frame( m1 = as.numeric( k_res1$pointwise),  m2 = as.numeric( k_res2$pointwise)) %>% 
  gather( model, value, m1:m2) %>% 
  ggplot( aes( x = value, fill = model, y = model )) + 
  ggridges::geom_density_ridges() 

  
model_combos[11, ]

model_combos <- read_rds('output/model_combos.rds')
load('data/temp_data/prepped_dfs.rda')
total <- length(model_combos)

base_model <- '~ 1'

# Model Parameters 
formX = as.formula(paste0 ('~ C')) ### Fixed effects design matrix (include climate as "C")
formZ = as.formula(paste0 ('~ 1'))

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
  split(., .$id)

model_combos <- model_combos[[1]]

id <- model_combos$id                      
sp <- model_combos$species
ad <- model_combos$ad
model <- model_combos$model
formC <- as.formula( paste0 ( '~-1 + ', model  ))  ### Climate effects design matrix 
fold <- model_combos$fold
formC

dat <- prepped_dfs[[sp]]  
hold <- unique( dat$yid[ dat$folds == fold ] )

model_combos
model_combos %>% read_csv('data/climate_window_descriptions.csv')

pblapply(model_combos,
         get_model_score,
         vr,
         model_file,
         prepped_rec_dfs,
         formZ = formZ,   
         formX = formX, 
         n_iter = n_iter,
         nthin = nthin)

