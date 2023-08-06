rm(list = ls() )
library(tidyverse)
source('code/analysis/functions.R')

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
sp <- "ARTR"

growth_cv <- read_csv('output/growth_models/cross_validation_growth_models.csv')

growth_delta_cv_stats <- growth_cv %>% 
  pivot_longer(AICc_clim:RMSE_null) %>% 
  separate( name, c('stat', 'model'), sep = '_') %>%  
  group_by( species, size_class, stat ) %>% 
  mutate( delta = value - value[model == 'null']) %>% 
  ungroup() %>% 
  filter( model != 'null') %>% 
  select( species, size_class, stat, model, delta ) %>% 
  unite( col = 'new',  c('stat', 'model'), sep = '_') %>% 
  pivot_wider(names_from = new, values_from = delta ) %>% 
  arrange( size_class, species ) %>% 
  split( .$size_class)
  

growth_delta_cv_stats$large  
growth_delta_cv_stats$small



# make a table out of this 
growth_cv %>% 
  arrange( size_class, species ) %>% View 


# rank models 
growth_model_ranks <- growth_cv %>% 
  pivot_longer(AICc_clim:RMSE_null) %>% 
  separate( name, c('stat', 'model'), sep = '_' ) %>% 
  group_by( size_class, species, stat) %>% 
  arrange(size_class, species, stat, value  ) %>% 
  filter( ! (size_class == 'small' & model == 'clim no intxn') ) %>% 
  mutate( model_rank = ifelse( stat == 'R2', rank(rev(value)), rank(value) )) %>% 
  arrange( size_class, species, stat, model_rank ) %>% 
  ungroup() 

growth_model_ranks %>%   
  filter(!stat %in% c('AICc', 'BIC', 'R2')) %>% 
  group_by( size_class, species, model ) %>%
  summarise( overall_model_rank = sum(model_rank ) ) %>% 
  arrange( size_class, species, overall_model_rank) %>% View 

growth_model_ranks %>%   
  filter(stat %in% c('AICc', 'BIC')) %>% 
  group_by( size_class, species, model ) %>%
  summarise( overall_model_rank = sum(model_rank ) ) %>% 
  arrange( size_class, species, overall_model_rank) %>% View 

# Survival models------------- #
survival_cv <- read_csv('output/survival_models/cross_validation_survival_models.csv')

surv_delta_cv_stats <- survival_cv %>% 
  pivot_longer(AUC_clim:DIC_null) %>% 
  separate( name, c('stat', 'model'), sep = '_') %>%  
  group_by( species, stat ) %>% 
  mutate( delta = value - value[model == 'null']) %>% 
  ungroup() %>% 
  filter( model != 'null') %>% 
  select( species, stat, model, delta ) %>% 
  unite( col = 'new',  c('stat', 'model'), sep = '_') %>% 
  pivot_wider(names_from = new, values_from = delta ) %>% 
  arrange( species ) 


surv_delta_cv_stats
# make a table out of this 
survival_cv 

# rank models 
surv_model_ranks <- survival_cv %>% 
  pivot_longer(AUC_clim:DIC_null) %>% 
  separate( name, c('stat', 'model'), sep = '_' ) %>% 
  group_by( species, stat) %>% 
  arrange(species, stat, value  ) %>% 
  mutate( model_rank = ifelse( stat == 'AUC', rank(rev(value)), rank(value) ) ) %>% 
  arrange( species, stat, model_rank ) %>% 
  ungroup() 

surv_model_ranks

surv_model_ranks %>%   
  filter(!stat %in% c('AICc', 'BIC', 'DIC')) %>% 
  group_by( species, model ) %>%
  summarise( overall_model_rank = sum(model_rank ) ) %>% 
  arrange( species, overall_model_rank) %>% View 

surv_model_ranks %>%   
  filter(stat %in% c('AICc', 'BIC', 'DIC')) %>% 
  group_by(  species, model ) %>%
  summarise( overall_model_rank = sum(model_rank ) ) %>% 
  arrange(  species, overall_model_rank) %>% View 

surv_model_ranks %>% 
  filter( stat == 'AUC') %>% 
  group_by( species, model) %>% 
  arrange( species, model_rank)
