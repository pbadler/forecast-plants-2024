rm(list = ls() )
library(tidyverse)
source('code/analysis/functions.R')
# For each species this script fits a: 
#  1. Survival model with size x climate interaction. 
#  2. Survival model without size x climate interaction 
#  3. Null survival model without climate 
#
#  Finally it saves training and testing data for each model.
# ----------------------------------------------------------  

split_year <- 2010 # Training testing split year 
size_cutoff <- -1  # log scale size cutoff between large and small 
quad_info <- read_csv( file = 'data/quad_info.csv')
daily_weather <- read_csv('data/temp/daily_weather_for_models.csv')

survival_windows <- read_csv('output/survival_models/top_survival_windows_by_deltaLogLoss.csv')
species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

sp <- 'ARTR'
for( sp in species_list){ 
  
  temp <- 
    survival_windows %>% 
    filter( species == sp , type == 'survival')
  
  surv <- prep_survival_for_climWin(species = sp, 
                                    last_year = 2020, # get all data 
                                    quad_info = quad_info) %>% 
    ungroup() %>% 
    mutate(size_class = ifelse( area0 > size_cutoff , "large", "small"))
  
  first_var <- temp[ temp$fit == 1, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')] 
  second_var <- temp[ temp$fit == 2, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')]
  
  clim1 <- getWindowAvg(daily_weather, 
                        var = str_remove( first_var$climate[1], '_scaled'), 
                        open = first_var$WindowOpen, 
                        close = first_var$WindowClose)
  
  clim2 <- getWindowAvg(daily_weather, 
                        var = str_remove( second_var$climate[1], '_scaled'), 
                        open = second_var$WindowOpen, 
                        close = second_var$WindowClose)
  
  temp_clim <- 
    clim1 %>% 
    left_join(clim2, by = c('Treatment', 'year')) %>% 
    mutate( !!first_var$climate[1]:= as.numeric( scale(.data[[str_remove(first_var$climate[1], '_scaled')]] ))) %>% 
    mutate( !!second_var$climate[1]:=as.numeric( scale(.data[[str_remove(second_var$climate[1], '_scaled')]]))) %>% 
    select( Treatment, year, 
            .data[[str_remove( first_var$climate[1], '_scaled')]], 
            .data[[str_remove( second_var$climate[1], '_scaled')]], 
            .data[[first_var$climate[1]]], 
            .data[[second_var$climate[1]]]) 

  temp_dat <- 
    surv %>% 
    filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>% 
    left_join(temp_clim, by = c('year', 'Treatment')) %>% 
    mutate( Split = ifelse( year <= split_year, 'Training', 'Testing')) %>% 
    split( .$Split)
  
  training <- 
    temp_dat$Training
  
  # Save testing and training data 
  temp_name <- paste0( sp, '_survival')
  write_csv( temp_dat$Training, paste0( 'data/temp/', temp_name, '_training_data.csv'))
  write_csv( temp_dat$Testing, paste0( 'data/temp/', temp_name, '_testing_data.csv'))
  
  # Save climate model 
  W.intra <- paste0( 'W.', sp)
  frm <- paste0( 'survives ~ area0*', first_var$climate[1], ' + area0*', second_var$climate[1], " + " , W.intra,   ' + (1|year/Group)' )
  my_mod <- glmer( formula = frm, 
                   data = training, 
                   family = 'binomial', 
                   control = glmerControl(optimizer = 'bobyqa'))
  
  saveRDS(my_mod, paste0( "output/survival_models/", temp_name, ".rds"))
  
  # Save null model 
  frm_null <- as.formula( paste0( 'survives ~ area0 + ', W.intra,   ' + (1|year/Group)' ))
  my_mod_null <- glmer( formula = frm_null, 
                   data = training, 
                   family = 'binomial', 
                   control = glmerControl(optimizer = 'bobyqa'))
  
  model_name <- paste0( sp, '_survival_null')
  saveRDS(my_mod_null, paste0( "output/survival_models/", temp_name, "_null.rds"))
 
  
  # -- Survival no intxn models -------------- # 
  temp <- 
    survival_windows %>% 
    filter( species == sp , type == 'survival_no_intxn')
  
  surv <- prep_survival_for_climWin(species = sp, 
                                    last_year = 2020, # get all data 
                                    quad_info = quad_info) %>% 
    ungroup() %>% 
    mutate(size_class = ifelse( area0 > size_cutoff , "large", "small"))
  
  first_var <- temp[ temp$fit == 1, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')] 
  second_var <- temp[ temp$fit == 2, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')]
  
  clim1 <- getWindowAvg(daily_weather, 
                        var = str_remove( first_var$climate[1], '_scaled'), 
                        open = first_var$WindowOpen, 
                        close = first_var$WindowClose)
  
  clim2 <- getWindowAvg(daily_weather, 
                        var = str_remove( second_var$climate[1], '_scaled'), 
                        open = second_var$WindowOpen, 
                        close = second_var$WindowClose)
  
  temp_clim <- 
    clim1 %>% 
    left_join(clim2, by = c('Treatment', 'year')) %>% 
    mutate( !!first_var$climate[1]:= as.numeric( scale(.data[[str_remove(first_var$climate[1], '_scaled')]] ))) %>% 
    mutate( !!second_var$climate[1]:=as.numeric( scale(.data[[str_remove(second_var$climate[1], '_scaled')]]))) %>% 
    select( Treatment, year, 
            .data[[str_remove( first_var$climate[1], '_scaled')]], 
            .data[[str_remove( second_var$climate[1], '_scaled')]], 
            .data[[first_var$climate[1]]], 
            .data[[second_var$climate[1]]]) 
  
  temp_dat <- 
    surv %>% 
    filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>% 
    left_join(temp_clim, by = c('year', 'Treatment')) %>% 
    mutate( Split = ifelse( year <= split_year, 'Training', 'Testing')) %>% 
    split( .$Split)
  
  training_no_intxn <- 
    temp_dat$Training
  
  # Save testing and training data 
  temp_name <- paste0( sp, '_survival_no_intxn')
  write_csv( temp_dat$Training, paste0( 'data/temp/', temp_name, '_training_data.csv'))
  write_csv( temp_dat$Testing, paste0( 'data/temp/', temp_name, '_testing_data.csv'))
  
  # Save climate model 
  W.intra <- paste0( 'W.', sp)
  frm <- paste0( 'survives ~ area0  + ', first_var$climate[1], ' + ', second_var$climate[1], " + " , W.intra,   ' + (1|year/Group)' )
  my_mod <- glmer( formula = frm, 
                   data = training, 
                   family = 'binomial', 
                   control = glmerControl(optimizer = 'bobyqa'))
  
  saveRDS(my_mod, paste0( "output/survival_models/", temp_name, ".rds"))
  
  
}
