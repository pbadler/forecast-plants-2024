rm(list = ls())
library(rstan)
library(tidyverse)

source('code/analysis/stan_data_functions.R')

model_table <- read_csv('output/IBM_model_table.csv')
quad_info <- read_csv('data/quad_info.csv')

load('code/figure_scripts/my_plotting_theme.Rdata')

compare_models <- 
  model_table %>% 
  gather( type, file, growth_data:survival_model) %>% 
  separate( type, c('vr', 'type')) %>% 
  spread( type, file )

compare_models$oos_lpd <- NA
compare_models$oos_lpd_irrigation <- NA
compare_models$oos_lpd_drought <- NA
compare_models$oos_lpd_control <- NA

i <- 1
for( i in 1:nrow(compare_models) ){   
  
  temp_fit <- readRDS(compare_models$model[i])
  
  oos_lpd <- sum( get_lpd(temp_fit, 'hold_log_lik') )
  
  hold_ll <- as.matrix(temp_fit, 'hold_log_lik')
  
  sum( log( colMeans( exp( hold_ll ) ) ) ) # all lpd 
  
  temp_dat <- readRDS(compare_models$data[i])
  
  Irrigation <- which( temp_dat$hold_quad_name %in% quad_info$QuadName[quad_info$Treatment == "Irrigation"] )
  Drought <- which( temp_dat$hold_quad_name %in% quad_info$QuadName[quad_info$Treatment == "Drought"] )
  Control <- which( temp_dat$hold_quad_name %in% quad_info$QuadName[quad_info$Treatment == "Control"] )
  
  Irrigation_lpd <- sum( log( colMeans( exp( hold_ll[, Irrigation] ) ) ) )
  Drought_lpd <- sum( log( colMeans( exp( hold_ll[, Drought] ) ) ) )
  Control_lpd <- sum( log( colMeans( exp( hold_ll[, Control] ) ) ) )
  
  compare_models$oos_lpd[i] <- oos_lpd 
  compare_models$oos_lpd_irrigation[i] <- Irrigation_lpd
  compare_models$oos_lpd_drought[i] <- Drought_lpd
  compare_models$oos_lpd_control[i] <- Control_lpd 

}

simple_table <- 
  compare_models %>% 
  gather( type, oos_lpd, starts_with('oos')) %>%
  select( spp, vr, climate, type, oos_lpd) %>% 
  mutate( climate = ifelse( climate, 'climate_model', 'null_model'))


simple_table %>%
  mutate( oos_lpd_label = str_extract(type, '[a-z]+$')) %>% 
  mutate( oos_lpd_label = ifelse( oos_lpd_label == 'lpd', 'total', oos_lpd_label)) %>% 
  ggplot( aes(x = oos_lpd_label, y = oos_lpd, color = climate)) + 
  geom_point() + 
  facet_wrap( vr ~ spp, scale = 'free', nrow = 3, ncol = 4) + 
  my_theme

