rm(list = ls() )
library(tidyverse)
source('code/analysis/functions.R')

split_year <- 2010 # Training testing split year 
size_cutoff <- -1  # log scale size cutoff between large and small 
quad_info <- read_csv( file = 'data/quad_info.csv')
daily_weather <- read_csv('data/temp/daily_weather_for_models.csv')

growth_windows <- read_csv('output/growth_models/top_growth_windows.csv')
sp <- 'ARTR'

temp <- 
  growth_windows %>% 
  filter( species == sp )

first_var <- temp[ temp$fit == 1, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')] 
second_var <- temp[ temp$fit == 2, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')]

clim1 <- getWindowAvg(daily_weather, var = first_var$climate[1], open = first_var$WindowOpen[1], close = first_var$WindowClose[1])
clim2 <- getWindowAvg(daily_weather, var = second_var$climate[1], open = second_var$WindowOpen[1], close = second_var$WindowClose[1])

temp_clim <- clim1 %>% 
  left_join(clim2, by = c('Treatment', 'year')) %>% 
  select( Treatment, year, first_var$climate[1], second_var$climate[1]) %>% 
  filter( Treatment == 'Control')

size <- prep_growth_for_climWin(species = 'ARTR', 
                                last_year = 2020, 
                                size_cutoff = -Inf, 
                                quad_info = quad_info)

TempDat <- size %>% 
  ungroup() %>%
  mutate(size_class = ifelse( area0 > size_cutoff , "large", "small")) %>% 
  ungroup() %>%
  left_join(temp_clim, by = c('year', 'Treatment')) %>% 
  mutate( Split = ifelse( year < split_year, 'Training', 'Testing')) %>% 
  split( .$Split)

training <- 
  TempDat$Training %>% 
  filter( size_class == 'large')

W.intra <- paste0( 'W.', sp)

frm <- paste0( 'area ~ area0*', first_var$climate[1], ' + area0*', second_var$climate[1], " + " , W.intra,   ' + (1|year/Group)' )
my_mod <- lmer( formula = frm, data = training, control = control_lmer, REML = F)

load('output/growth_models/ARTR_growth_mer_monthly_ClimWin.rda')
ARTR_growth_monthly_ClimWin$ClimWinFit2$combos
winMod <- ARTR_growth_monthly_ClimWin$ClimWinFit2[[1]]$BestModel

winMod@frame %>% nrow 
training %>% select(area, area0, climate  )

my_mod@frame %>% nrow 

fixef(winMod)
fixef(my_mod)

load( 'output/growth_models/ARTR_growth_mer_monthly_ClimWin.rda')

fixef( ARTR_growth_monthly_ClimWin$ClimWinFit2[[1]]$BestModel )
fixef( read_rds('output/growth_models/ARTR_growth.rds'))
fixef( my_mod)


