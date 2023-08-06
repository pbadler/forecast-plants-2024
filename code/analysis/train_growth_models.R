rm(list = ls() )
library(tidyverse)
library(lubridate)
source('code/analysis/functions.R')
# For each species this script fits a: 
#  1. Growth model with size x climate interaction. 
#  2. Growth model without size x climate interaction 
#  3. Null growth model without climate 
#  4. Growth model for small plants 
#  5. Null growth model for small plants 
#  
#  Finally it saves training and testing data for each model.
# ----------------------------------------------------------  

split_year <- 2010 # Training testing split year 
size_cutoff <- -1  # log scale size cutoff between large and small 
quad_info <- read_csv( file = 'data/quad_info.csv')
daily_weather <- read_csv('data/temp/daily_weather_for_models.csv')
growth_windows <- read_csv('output/growth_models/top_growth_windows_by_deltaMSE.csv')

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

for( sp in species_list){ 
  
  temp <- 
    growth_windows %>% 
    filter( species == sp, 
            type == 'growth')
  
  
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

  temp_clim <- clim1 %>% 
    left_join(clim2, by = c('Treatment', 'year')) %>% 
    mutate( !!first_var$climate[1]:= as.numeric( scale(.data[[str_remove(first_var$climate[1], '_scaled')]] ))) %>% 
    mutate( !!second_var$climate[1]:=as.numeric( scale(.data[[str_remove(second_var$climate[1], '_scaled')]]))) %>% 
    dplyr::select( Treatment, year, 
            all_of(str_remove( first_var$climate[1], '_scaled')), 
            all_of(str_remove( second_var$climate[1], '_scaled')), 
            all_of(first_var$climate[1]), 
            all_of(second_var$climate[1]))
  
  size <- prep_growth_for_climWin(species = sp, 
                                  last_year = 2020, # get all data 
                                  size_cutoff = -Inf, # get all data 
                                  quad_info = quad_info) %>% 
    mutate(size_class = ifelse( area0 > size_cutoff , "large", "small"))
  
  # Split training and testing data ------------- 
  temp_dat <- 
    size %>% 
    filter( size_class == 'large') %>% 
    filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>% 
    ungroup() %>%
    left_join(temp_clim, by = c('year', 'Treatment')) %>% 
    mutate( Split = ifelse( year <= split_year, 'Training', 'Testing')) %>% 
    split( .$Split)
  # ---------------------------------------------
  
  temp_name <- paste0( sp, '_growth')
  
  # Save Training and Testing data 
  stopifnot( max( temp_dat$Training$year ) <= split_year )
  stopifnot( min( temp_dat$Testing$year ) > split_year )
  
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
  
  # no intxn models 
  temp <- 
    growth_windows %>% 
    filter( species == sp, 
            type == 'growth_no_intxn')
  
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
  
  
  temp_clim <- clim1 %>% 
    left_join(clim2, by = c('Treatment', 'year')) %>% 
    mutate( !!first_var$climate[1]:= as.numeric( scale(.data[[str_remove(first_var$climate[1], '_scaled')]] ))) %>% 
    mutate( !!second_var$climate[1]:=as.numeric( scale(.data[[str_remove(second_var$climate[1], '_scaled')]]))) %>% 
    dplyr::select( Treatment, year, 
                 all_of(str_remove( first_var$climate[1], '_scaled')), 
                 all_of(str_remove( second_var$climate[1], '_scaled')), 
                 all_of(first_var$climate[1]), 
                 all_of(second_var$climate[1]))
  
  
  # split training and testing ------------- # 
  temp_dat <- 
    size %>% 
    filter( size_class == 'large') %>% 
    filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>% 
    ungroup() %>%
    left_join(temp_clim, by = c('year', 'Treatment')) %>% 
    mutate( Split = ifelse( year <= split_year, 'Training', 'Testing')) %>% 
    split( .$Split)  
  
  temp_name <- paste0( sp, '_growth_no_intxn')
  
  write_csv( temp_dat$Training, paste0( 'data/temp/', temp_name, '_training_data.csv'))
  write_csv( temp_dat$Testing, paste0( 'data/temp/', temp_name, '_testing_data.csv'))
  
  # --------------------------------------------------------- #
  training_no_intxn <- temp_dat$Training %>% filter( size_class == 'large')
  
  # Climate informed model 
  frm <- as.formula( paste0( 'area ~ ', 'area0  +' , first_var$climate[1], ' + ', second_var$climate[1], " + " , W.intra,   ' + (1|year/Group)' ))
  
  my_mod <- lmer( formula = frm, data = training_no_intxn, control = control_lmer, REML = T)
  saveRDS(my_mod, paste0( "output/growth_models/", temp_name, ".rds"))
  
  # small plant models 
  temp <- 
    growth_windows %>% 
    filter( species == sp, 
            type == 'growth_small')
  
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
  
  temp_clim <- clim1 %>% 
    left_join(clim2, by = c('Treatment', 'year')) %>% 
    mutate( !!first_var$climate[1]:= as.numeric( scale(.data[[str_remove(first_var$climate[1], '_scaled')]] ))) %>% 
    mutate( !!second_var$climate[1]:=as.numeric( scale(.data[[str_remove(second_var$climate[1], '_scaled')]]))) %>% 
    dplyr::select( Treatment, year, 
                 all_of(str_remove( first_var$climate[1], '_scaled')), 
                 all_of(str_remove( second_var$climate[1], '_scaled')), 
                 all_of(first_var$climate[1]), 
                 all_of(second_var$climate[1]))
  
  
  # split training and testing ------------- # 
  temp_dat <- 
    size %>% 
    filter( size_class == 'small') %>% 
    filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>% 
    ungroup() %>%
    left_join(temp_clim, by = c('year', 'Treatment')) %>% 
    mutate( Split = ifelse( year <= split_year, 'Training', 'Testing')) %>% 
    split( .$Split)  
  
  temp_name <- paste0( sp, '_growth_small')
  
  write_csv( temp_dat$Training, paste0( 'data/temp/', temp_name, '_training_data.csv'))
  write_csv( temp_dat$Testing, paste0( 'data/temp/', temp_name, '_testing_data.csv'))
  
  # --------------------------------------------------------- #
  training_small <- temp_dat$Training %>% filter( size_class == 'small')
  
  # Climate informed model 
  frm <- as.formula( paste0( 'area ~ ', first_var$climate[1], ' + ', second_var$climate[1], " + " , W.intra,   ' + (1|year)' ))
  my_mod <- lmer( formula = frm, data = training_small, control = control_lmer, REML = T)
  saveRDS(my_mod, paste0( "output/growth_models/", temp_name, ".rds"))
  
  # Null model  
  frm_null <- as.formula(  frm <- as.formula( paste0( 'area ~ ', W.intra,   ' + (1|year)' )))
  my_mod_null <- lmer( formula = frm_null, data = training_small, control = control_lmer, REML = T)
  
  saveRDS(my_mod_null, paste0( "output/growth_models/", temp_name, "_null.rds"))
  
    
}
