par(mfrow= c(1,1))
library(tidyverse)
library(imputeTS)
dat <- read_csv('lib/driversdata/data/idaho/speciesData/ARTR/ARTR_genet_xy.csv')

dat <- 
  dat %>% 
  mutate(  pid = paste( quad, trackID, sep = '_'))

mnyear <- min(dat$year)
mxyear <- max(dat$year)
quad <- unique( dat$quad)
pid <- unique( dat$pid )
year <- mnyear:mxyear

dat <- 
  expand.grid( year = mnyear:mxyear, pid = unique(pid)) %>% 
  left_join(dat, by = c('pid', 'year'))  %>% 
  arrange(pid, year) 

my_dat <- 
  dat %>% 
  group_by( pid ) %>% 
  mutate( age = ifelse( year == min(year[!is.na(area)]), 1, age )) %>% 
  dplyr::select( pid, age, year, area )


my_dat <- my_dat %>% mutate( n = sum(!is.na(area))) %>% arrange(desc(n))

my_ts <- split(my_dat, my_dat$pid) %>% lapply( function(x) ts(log( x$area)) )
my_ts <- my_ts[  order( unlist( lapply( my_ts , function(x) sum(!is.na(x))) ) , decreasing = T) ]

test <- my_ts[[6]]
test_interp <- na_interpolation(test, 'linear')
missing <- which(is.na(test))
ip <- test_interp[ which(is.na(test))] 
plot( test, type = 'n')
points(test_interp, type = 'l')
points( na_interpolation(test,  'stine'), type = 'l', col = 'green')
points(test, type = 'b')
points( missing, ip , pch = 20)

has_values <- unlist( lapply( my_ts , function(x)sum(!is.na(x)) > 1 )  )

test_intpd <- lapply( my_ts[has_values], na_interpolation, 'linear')

max(my_dat$year)
min(my_dat$year)

plot(test_intpd[[1]])
test <- arima(test_intpd[[1]][0:36], c(1, 1, 2), seasonal = c(0, 0, 0))
#plot(resid(test), type = 'n')
points(resid(test), type = 'p')

# find short sequences of missing values
# interpolate these
# fit time series models to each plant with long sequence
# look at residuals 
# look for year effects in the residuals 