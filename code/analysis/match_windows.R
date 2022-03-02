library(lubridate)
library(brms)
library(tidyverse)

species <- 'HECO'

load( 'data/temp_data/HECO_growth_weekly_ClimWin.rda')

model_index <- 2  # This is the model with log VWC 

test_data <- HECO_growth_weekly_ClimWin[[model_index]]$BestModelData

select_window <- HECO_growth_weekly_ClimWin$combos[model_index, ]

daily_vwc <- 
  readRDS( 'data/temp_data/daily_swVWC_treatments.RDS') %>% 
  filter( Treatment == 'Control' ) %>% 
  arrange( date )

# See if I can recreate test data vwc
quad_info <- read_csv('data/quad_info.csv')
last_year <- 2010
intra <- paste0( 'W.' , species )
dat <- readRDS( paste0 ( 'data/temp_data/', paste0( species, '_size.RDS')))
dat <- dat[ dat$year < last_year, ]

mnyear <- min(dat$year)
mxyear <- max(dat$year)
pid <- unique( dat$pid )

dat <- 
  expand.grid( year = mnyear:mxyear, pid = unique(pid)) %>% 
  left_join(dat, by = c('pid', 'year'))  %>% 
  arrange(pid, year) 

my_dat <- 
  dat %>% 
  dplyr::select(quad, pid, age, year, area, starts_with('W.'))  %>% 
  mutate( W.total = W.ARTR + W.HECO + W.POSE + W.PSSP - eval(parse(text = intra)), 
          W.intra = eval(parse(text = intra)))

growth <- 
  my_dat %>% 
  group_by( pid ) %>% 
  arrange(pid, year) %>% 
  mutate( area = log(area)) %>% 
  mutate( area0 = lag(area)) %>% 
  mutate( growth = area - area0 ) %>% 
  filter( !is.na(growth),  
          is.finite(growth), 
          area > -1.38,                   # filter to larger plants as more reliable indicators of growth
          area0 > -1.38, year > 1925) %>%
  mutate( date = paste0( '15/06/', year)) # date for climwin 

growth <- 
  growth %>% 
  left_join(quad_info, by = 'quad') %>%
  mutate( group_year = paste( Group , year, sep = '_'))

bdate <- dmy( growth$date  )
cdate <- ymd(daily_vwc$date)

cweek <- week( cdate)
cweek[ which(cweek==53) ] <- 52 # Set extra few days per year in week 53 to week 52
bweek <- week(bdate)

cyear <- year( cdate ) - min(year(cdate))
cweek_no <- cweek + 52*cyear

byear <- year( bdate ) - min(year(cdate))
bweek_no <- bweek + 52*byear

growth$bweek_no <- bweek_no
window_range <- (select_window$WindowClose:select_window$WindowOpen)


daily_vwc$cweek_no <- cweek_no

new_VWC <- NA
u_bweek_no <- unique(bweek_no)

for( i in 1:length(u_bweek_no)){ 

  new_VWC[i] <- log( mean( daily_vwc$VWC[ (u_bweek_no[i] - cweek_no) %in% window_range ] ))
    
}

new_growth <- 
  growth %>%
  left_join(
    data.frame( bweek_no = u_bweek_no , climate = new_VWC) ) 

all.equal(test_data$`log(climate)`, new_growth$climate)

plot( new_growth$climate, test_data$`log(climate)`) 



