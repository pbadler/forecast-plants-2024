rm(list = ls())
library(future)
library(loo)
library(brms)
library(tidyverse)
library(forecast)

source('code/analysis/stan_data_functions.R')

rmse <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  sqrt(mean((yrep_mean - y)^2))
}

testing <- F
n_iter <- 5000 

small <- -0.5
lc <- -0.5
k <- 10
species <- 'PSSP'
vr <- 'growth'
dat <- readRDS('data/temp_data/ARTR_growth_survival_dataframe.RDS')

dat <- 
  dat %>% 
  select( - W ) %>% 
  filter( Period == 'Historical') %>% 
  filter( !is.na(logarea.t1 ), !is.na(logarea.t0)) 

if(testing) { 
  dat <-
    dat %>% 
    group_by(yid) %>% 
    sample_frac(0.25, replace = F) %>% 
    ungroup() 
}

dat$Y <- dat$logarea.t1
dat$size <- dat$logarea.t0
dat$small <- as.numeric(dat$logarea.t0 < small)

dat <- 
  dat %>% 
  mutate( W.total = rowSums( select( . , starts_with('W')))) %>% 
  mutate( W.intra = rowSums( select( . , !!(paste0('W.', species))))) %>% 
  mutate( W.inter = W.total - W.intra)

dat$censored <- 0
dat$censored <- (dat$Y <= lc)*-1 
dat$plant_id <- paste( dat$quad, dat$trackID, sep = '_') 

plot( dat$Y )
abline( h = lc )
plot(dat$size)
abline( h = small )

plot(dat$size, dat$Y)


model_combos <- read_rds('output/model_combos.rds') %>% distinct(model)
C_spec <- paste(model_combos$model[2:8], collapse = '+' )
group_spec <- '(1|yid) + (1|plant_id)'
LHS <- 'Y | cens(censored)'
RHS <- paste( 'size + W.inter + W.intra', C_spec, group_spec, sep = ' + ' )
my_form <- paste(LHS, RHS, sep = ' ~ ')

scaled_covs <- 
  dat %>% 
  select( Y, size, starts_with('C.T'), starts_with('C.V'), starts_with('W.'))

dat <- 
  dat %>% 
  select( year, age, quad, trackID, censored ) %>% 
  bind_cols(
    scaled_covs %>% 
      select( Y, size, starts_with('C.T'), starts_with('C.V'), starts_with('W')) %>% 
      mutate_all(scale))


testdat <- 
  dat %>% 
  filter( censored == 0  ) %>% 
  mutate( plantID = paste( quad, trackID, sep = '_')) %>% 
  select( - starts_with('C'), -starts_with('W')) 

pts <- expand.grid( plantID = unique(testdat$plantID) , year = (min(testdat$year) - 5):max(testdat$year + 5) ) %>% 
  left_join(testdat) %>% 
  arrange( plantID, year) %>% 
  group_by( plantID ) %>% 
  mutate( Y = scale(Y))

pts %>% 
  ggplot( aes( x = year, y = Y, group = plantID))  + geom_line() + 
  geom_hline(aes( yintercept = 0))

myts <- ts( as.numeric( pts$Y ) )

plot(myts)
length(myts)

library(forecast)
arfit <- auto.arima(myts,seasonal = F, allowdrift = F)

plot(arfit)
plot( residuals(arfit), type ='l')
summary(arfit)
