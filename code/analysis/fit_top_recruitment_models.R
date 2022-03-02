rm(list = ls())

library(rstan)
library(tidyverse)
library(loo)

source('code/analysis/stan_data_functions.R')

vr <- 'recruitment'
stan_model_file <- 'code/analysis/recruitment.stan'

testing <- F
if( testing ){ 
  
  # STAN pars -------------- 
  ncores <- 1 
  niter <- 1000
  nchains <- 1 
  nthin <- 5
  
}else{
  
  # STAN pars -------------- 
  ncores <- 4 
  niter <- 2000
  nchains <- 4 
  nthin <- 5
  
}

formX = as.formula('~ C') ### Fixed effects design matrix (include climate as "C")
# ------------------------------------------

# set up climate variable table --------------------------# 
stan_mods <- read_csv('output/model_ranks.csv')

top_mods <- 
  stan_mods %>% 
  group_by( spp, vr ) %>% 
  filter(oos_lppd == max(oos_lppd), vr == vr) %>% 
  select( vr, spp, climate_window)

model_list <- expand.grid( 
  spp = unique( top_mods$spp) , 
  vr = vr, 
  model = c('top_model', 'none'))

model_list <- 
  model_list %>% 
  left_join(top_mods) %>% 
  mutate( climate_window = ifelse(model == 'none', 'none', climate_window)) 

# --------------------------------------------------------- #
nsp <- length(unique( model_list$spp))

model_list$adapt_delta <- c(0.98, 0.98, 0.8, 0.8)[nsp]
model_list$formX <- list( formX  )

formXNULL <- update(formX,  ~ . - C)

model_list <- 
  model_list %>% 
  mutate(formX = ifelse( climate_window == "none", list( formXNULL), formX  )) %>% 
  filter(!is.na(formX)) %>% 
  distinct( spp, vr, adapt_delta, climate_window, formX )

for(i in 1:nrow(model_list)){ 
  
  # choose species 
  sp <- model_list$spp[i]
  ad <- model_list$adapt_delta[i]
  fx <- model_list$formX[[i]]
  window <- model_list$climate_window[i]
  
  dat_file <- paste0('data/temp_data/', sp, '_recruitment_dataframe.RDS')
  dat <- readRDS(dat_file)
  
  dat <- 
    dat %>% 
    ungroup %>% 
    rowwise( ) %>% 
    mutate( total_basal_cover = cov.HECO + cov.POSE + cov.PSSP) %>% 
    mutate( total_open = 100*100 - total_basal_cover) %>% 
    mutate( l_open = log(total_open))
  
  # dat$P1 <- dat[ , paste0('cov.', sp) ]
  # dat$P2 <- dat[ , paste0('Gcov.', sp)]
  # 
  # dat$P1_inter <- dat$all_cover - dat$P1
  # 
  # dat$P1_inter <- scale( sqrt( dat$P1_inter))
  # dat$P2 <- scale( sqrt( dat$P2 ))
  # dat$open <- log(10000 - dat$all_cover)
  # 
  # dat$open
  
  if ( window != 'none' ){   
    moist <- paste0( 'C.VWC.', window)
    therm <- paste0( 'C.T.', window )
    
    # get climate effects 
    formC <- as.formula( paste0 ( '~-1 + ', moist, '*', therm  ))  ### Climate effects design matrix 
  }else if( window == 'none'){ 
    
    formC <- as.formula( '~-1')  
  }
  
  hold <- unique( dat$yid [ dat$Period == 'Modern'] ) 
  
  dl <- process_recruitment_data(dat = dat, 
                                 formX = formX, 
                                 formC = formC, 
                                 center = T,
                                 hold = hold, 
                                 IBM = 1)
  
  mod <- rstan::stan_model(stan_model_file) # load stan model 
  
  print( paste( '### ---- species', sp, '; climate window', window, '--------------------##'))
  print( paste( '### ---- working on model', i, 'of', nrow(model_list),' -------------###' ))

  fit <- rstan::sampling(mod, 
                          data = dl, 
                          chains = nchains, 
                          iter = niter, 
                          cores = ncores,
                          thin = nthin,
                          pars = c('hold_log_lik', 'hold_SSE', 'beta', 'IBM_Y_hat'), 
                          control = list(adapt_delta = ad), 
                          refresh = -1)
  
  saveRDS(dl, file = paste0( 'output/stan_fits/', sp, '_', vr, '_', window, '_data.RDS'))
  saveRDS(fit, file = paste0( 'output/stan_fits/', sp, '_', vr, '_', window, '_model.RDS'))
  
}
