par(mfrow= c(1,1))
library(tidyverse)
library(nlraa)
library(fitdistrplus)

dat <- read_csv('lib/driversdata/data/idaho/speciesData/POSE/POSE_genet_xy.csv')

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


dat %>% 
  ggplot( aes( x = year, y = area, group = pid)) + 
  geom_line(data = dat[!is.na(dat$area), ], aes( x = year, y = area, group = pid), color = 'lightblue', alpha = 0.5) +
  geom_line(color = 'darkgray') + 
  geom_point( size = 2, shape = 1, alpha = 0.5) 




dat2 <- 
  dat %>% 
  group_by( pid ) %>% 
  mutate( age = ifelse( year == min(year[!is.na(area)]), 1, age ))


dat2 %>% 
  group_by( pid ) %>% 
  filter( area[ which(age == 1) ] < 0.5 ) %>% 
  ggplot( aes( x = age, y = area, group = pid)) + 
  geom_line(color = 'darkgray') + 
  geom_point( size = 2, shape = 1, alpha = 0.5) + 
  scale_y_log10() + 
  facet_wrap(~quad)

dat3 <- 
  dat2 %>% 
  group_by( pid ) %>% 
  filter( area[ which(age == 1) ] < 0.4 ) 


dat4 <- 
  dat3 %>% 
  group_by( pid ) %>% 
  mutate( with_data = sum( !is.na(area)) ) %>% 
  arrange( desc(with_data)) %>% 
  filter( with_data > 10) %>% 
  group_by( pid) %>%
  mutate( size = log(area)) %>% 
  mutate(std_size = size/(max(size, na.rm = T)))


dat4 %>% 
  ungroup() %>% 
  ggplot( aes( x = age, y = std_size)) + 
  geom_point() +
  geom_smooth(aes( group = pid), span = 3, se= F, size = 0.5, alpha = 0.4) + 
  facet_wrap( ~ with_data, scales = 'free_y')


my_fun <- function( x , fit ){ 
  x <- data.frame( age = x); 
  predict( fit, newdata = x ) 
}

i <- 1
fits <- list(NA)

library(brms)
prior1 <- prior(normal(0, 1), nlpar = "alpha") +
  prior(normal(1, 1), nlpar = "beta") + 
  prior(normal( 1, 1), nlpar = 'theta')

simple_dat <- 
  dat4 %>% 
  dplyr::select( year, pid, quad, age, std_size )  %>% 
  filter( !is.na(age), !is.na(std_size))


qfit <- brm(std_size ~ poly(age, 2) + (1|pid), data = simple_dat, family = 'normal') 
plot(qfit)
stanplot(qfit)
my_preds <- predict(qfit)

data.frame(simple_dat[!is.na(simple_dat$std_size), ], data.frame(my_preds)) %>% 
  ggplot( aes( x = age, y= std_size)) + 
  geom_point() + 
  geom_line(aes(x = age, y= Estimate, group = pid)) + 
  facet_wrap( ~with_data)





for( i in 1:length(unique(dat4$pid)) ){ 
  
  temp_pid <- unique(dat4$pid)[i]
  
  temp <- 
    dat4 %>% 
    filter( pid == temp_pid, !is.na(area)) %>% 
    arrange( age) 
  

  gfit <- try( nls(area ~ theta*dgamma( age, shape, rate),
              start = list( theta = 10, shape = 50, rate = 6 ), data = temp , 
              lower = c(0, 0, 0), 
              algorithm = 'port') )
  
  if( class( gfit ) != 'nls'){ 
    inits <- coef(fits[[i-1]][[1]])
    gfit <- try( nls(area ~ theta*dgamma( age, shape, rate),
                     start = inits, data = temp , 
                     lower = c(0, 0, 0), 
                     algorithm = 'port') )
  }
  
  
  nfit <- try( nls(area ~ theta*dnorm( age, center, scale),
              start = list( theta = 10, center = 10, scale = 10 ), data = temp , 
              lower = c(0, 0, 0), 
              algorithm = 'port') )
  
  if( class( nfit ) != 'nls'){ 
    inits <- coef(fits[[i-1]][[2]])
    gfit <- try( nls(area ~ theta*dnorm( age, center, scale),
                     start = inits, data = temp , 
                     lower = c(0, 0, 0), 
                     algorithm = 'port') )
  }
  
  tfit <- try( nls(area ~ theta*dt((age - center)/scale, spread), 
              start = list( theta = 10, center = 12, scale = 5 , spread = 0.5), data = temp , 
              lower = c(0.5, 0.5, 0.5, 1e-1), 
              algorithm = 'port') )
  
  if( class( tfit ) != 'nls'){ 
    inits <- coef(fits[[i-1]][[3]])
    gfit <- try( nls(area ~  theta*dt((age - center)/scale, spread),
                     start = inits, data = temp , 
                     lower = c(0, 0, 0), 
                     algorithm = 'port') )
  }
  fits[[i]] <- list(gfit, nfit, tfit, data = temp )
  
}

