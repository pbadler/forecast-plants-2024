rm(list = ls() )
library(tidyverse)
source('code/analysis/functions.R')

split_year <- 2010 # Training testing split year 
size_cutoff <- -1  # log scale size cutoff between large and small 
quad_info <- read_csv( file = 'data/quad_info.csv')
daily_weather <- read_csv('data/temp/daily_weather_for_models.csv')

survival_windows <- read_csv('output/survival_models/top_survival_windows.csv')
species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

sp <- 'ARTR'
for( sp in species_list){ 
  
  temp <- 
    survival_windows %>% 
    filter( species == sp )
  
  first_var <- temp[ temp$fit == 1, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')] 
  second_var <- temp[ temp$fit == 2, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')]
  
  clim1 <- getWindowAvg(daily_weather, var = first_var$climate[1], open = first_var$WindowOpen, close = first_var$WindowClose)
  clim2 <- getWindowAvg(daily_weather, var = second_var$climate[1], open = second_var$WindowOpen, close = second_var$WindowClose)
  
  #first_var$
  temp_clim <- clim1 %>% 
    left_join(clim2, by = c('Treatment', 'year')) %>% 
    select( Treatment, year, first_var$climate[1], second_var$climate[1]) %>% 
    mutate( !!paste0( first_var$climate[1],"_scaled"):=scale(.data[[first_var$climate[1]]])) %>% 
    mutate( !!paste0( second_var$climate[1], "_scaled"):=scale(.data[[second_var$climate[[1]]]]))
  

  surv <- prep_survival_for_climWin(species = sp, 
                                  last_year = 2020, # get all data 
                                  quad_info = quad_info) %>% 
    ungroup() %>% 
    mutate(size_class = ifelse( area0 > size_cutoff , "large", "small"))
  

  temp_dat <- 
    surv %>% 
    filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>% 
    left_join(temp_clim, by = c('year', 'Treatment')) %>% 
    mutate( Split = ifelse( year < split_year, 'Training', 'Testing')) %>% 
    split( .$Split)
  
  training <- 
    temp_dat$Training
  
  
  # Save testing and training data 
  temp_name <- paste0( sp, '_survival')
  temp_dat$Training <- temp_dat$Training %>% mutate_at( vars( contains('scaled')), as.numeric) 
  temp_dat$Testing <- temp_dat$Testing %>% mutate_at( vars( contains('scaled')), as.numeric) 
  write_csv( temp_dat$Training, paste0( 'data/temp/', temp_name, '_training_data.csv'))
  write_csv( temp_dat$Testing, paste0( 'data/temp/', temp_name, '_testing_data.csv'))
  
  # Save climate model 
  W.intra <- paste0( 'W.', sp)
  frm <- paste0( 'survives ~ area0*', first_var$climate[1], '_scaled + area0*', second_var$climate[1], "_scaled + " , W.intra,   ' + (1|year/Group)' )
  my_mod <- glmer( formula = frm, 
                   data = training, 
                   family = 'binomial', 
                   control = glmerControl(optimizer = 'bobyqa'))
  
  model_name <- paste0( sp, '_survival')
  saveRDS(my_mod, paste0( "output/survival_models/", model_name, ".rds"))
  
  # Save null model 
  frm_null <- as.formula( paste0( 'survives ~ area0 + ', W.intra,   ' + (1|year/Group)' ))
  my_mod_null <- glmer( formula = frm_null, 
                   data = training, 
                   family = 'binomial', 
                   control = glmerControl(optimizer = 'bobyqa'))
  
  model_name <- paste0( sp, '_survival_null')
  saveRDS(my_mod_null, paste0( "output/survival_models/", model_name, ".rds"))
  
}
