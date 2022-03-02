rm(list = ls())
library(tidyverse)

source('code/analysis/stan_data_functions.R')

# --------------------------
k <- 10
species <- c('ARTR', 'HECO', 'POSE', 'PSSP')

# Model Parameters 
small <- -1                 ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(paste0 ('~ size + small + W.intra + W.inter + C')) ### Fixed effects design matrix (include climate as "C")

left_cut <- data.frame( 
  species = species, 
  lc = c(0, -1, -1, -1), 
  ad = c(0.98, 0.98, 0.8, 0.8)
) # for censored data 

# set up model selection table --------------------------# 
load('data/temp_data/climate_combos.RData')

nvars <- length(T_combos)

climate_effects <- paste0( 'C.T.', 1:nvars, '*', 'C.VWC.', 1:nvars)
climate_effects <- c('NULL', climate_effects)
climate_effects <- climate_effects

model_combos <- 
  expand.grid(species = species, model = climate_effects, fold = 1:k) %>% 
  arrange( species, model, fold) %>% 
  mutate( id = row_number() ) %>% 
  select( id, species, model, fold ) %>% 
  left_join(left_cut)

gs_filenames <- dir(path = 'data/temp_data', pattern = '.*_growth_survival_dataframe.RDS', full.names = T)

prepped_dfs <- 
  mapply(x = as.character( left_cut$species), 
         y = left_cut$lc, 
         FUN = function(x, y, ... ) scale_and_fold(species = x, lc = y,...), 
         k = k, vr = 'growth', small = small, filenames = gs_filenames, SIMPLIFY = T)


rec_filenames <- dir(path = 'data/temp_data', pattern = '.*_recruitment_dataframe.RDS', full.names = T)

prepped_rec_dfs <- 
  mapply(x = as.character( left_cut$species), 
         y = left_cut$lc, 
         FUN = function(x, y, ... ) scale_and_fold(species = x, lc = y,...), 
         k = k, 
         vr = 'recruitment', 
         filenames = rec_filenames, SIMPLIFY = F)



save(prepped_dfs, prepped_rec_dfs,  file = 'data/temp_data/prepped_dfs.rda')

saveRDS( model_combos, file = 'output/model_combos.rds')
saveRDS( left_cut, file = 'output/left_cut.rds')

