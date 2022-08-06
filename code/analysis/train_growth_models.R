rm(list = ls() )
library(tidyverse)
source('code/analysis/functions.R')

split_year <- 2010 # Training testing split year 
size_cutoff <- -1  # log scale size cutoff between large and small 
quad_info <- read_csv( file = 'data/quad_info.csv')
daily_weather <- read_csv('data/temp/daily_weather_for_models.csv')

growth_windows <- read_csv('data/temp/top_growth_windows.csv')



species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
sp <- 'ARTR'
for( sp in species_list){ 
  
  temp <- 
    growth_windows %>% 
    filter( species == sp )
  
  first_var <- temp[ temp$fit == 1, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')] 
  second_var <- temp[ temp$fit == 2, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')]
  
  clim1 <- getWindowAvg(daily_weather, var = first_var$climate[1], open = first_var$WindowOpen, close = first_var$WindowClose)
  clim2 <- getWindowAvg(daily_weather, var = second_var$climate[1], open = second_var$WindowOpen, close = second_var$WindowClose)
  
  temp_clim <- clim1 %>% 
    left_join(clim2, by = c('Treatment', 'year')) %>% 
    select( Treatment, year, first_var$climate[1], second_var$climate[1])
  
  size <- prep_growth_for_climWin(species = sp, 
                                  last_year = 2020, # get all data 
                                  size_cutoff = -Inf, # get all data 
                                  quad_info = quad_info) %>% 
    mutate(size_class = ifelse( area0 > size_cutoff , "large", "small"))
  
  temp_dat <- 
    size %>% 
    filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>% 
    ungroup() %>%
    left_join(temp_clim, by = c('year', 'Treatment')) %>% 
    mutate( Split = ifelse( year < split_year, 'Training', 'Testing')) %>% 
    split( .$Split)
  
  temp_name <- paste0( sp, '_growth')
  
  # Save Training and Testing data 
  write_csv( temp_dat$Training, paste0( 'data/temp/', temp_name, '_training_data.csv'))
  write_csv( temp_dat$Testing, paste0( 'data/temp/', temp_name, '_testing_data.csv'))
  
  
  training_large <- 
    temp_dat$Training %>% 
    filter( size_class == 'large')
  
  W.intra <- paste0( 'W.', sp)
  
  # Climate informed model 
  frm <- as.formula( paste0( 'area ~ area0*', first_var$climate[1], ' + area0*', second_var$climate[1], " + " , W.intra,   ' + (1|year/Group)' ))
  my_mod <- lmer( formula = frm, data = training_large, control = control_lmer, REML = T)
  model_name <- paste0( sp, '_growth')
  saveRDS(my_mod, paste0( "output/growth_models/", model_name, ".rds"))

  # Null model  
  frm_null <- as.formula(  frm <- as.formula( paste0( 'area ~ area0 + ', W.intra,   ' + (1|year/Group)' )))
  my_mod_null <- lmer( formula = frm_null, data = training_large, control = control_lmer, REML = T)
  model_name <- paste0( sp, '_growth_null')
  saveRDS(my_mod_null, paste0( "output/growth_models/", model_name, ".rds"))
  
}
