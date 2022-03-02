# refit top growth models 

rm(list = ls())

library(rstan)
library(tidyverse)
library(loo)

source('code/analysis/stan_data_functions.R')

vr <- 'growth'
stan_model_file <- 'code/analysis/growth.stan'

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

small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(paste0 ('~ size + small + W.intra + W.inter + C')) ### Fixed effects design matrix (include climate as "C")
formE = as.formula(~ size)  ### For growth model, size dependent variance design matrix  
# ------------------------------------------

# set up climate variable table --------------------------# 
stan_mods <- read_csv('output/model_ranks.csv')

top_mods <- 
  stan_mods %>% 
  group_by( spp, vr ) %>% 
  filter(oos_lppd == max(oos_lppd), vr == vr) %>% 
  select( vr, spp, climate_window)

model_list <- expand.grid( spp = unique( top_mods$spp), 
                           vr = vr, 
                           model = c('top_model', 'none'))

model_list <- 
  model_list %>% 
  left_join(top_mods) %>% 
  mutate( climate_window = ifelse(model == 'none', 'none', climate_window)) 

# --------------------------------------------------------- #
nsp <- length(unique( model_list$spp))

model_list$adapt_delta <- c(0.9, 0.9, 0.8, 0.8)[nsp]
model_list$left_cut <- c(-1, -1.3, -1.3, -1.3)[nsp]
model_list$formX <- list( formX  )

formXNULL <- update(formX,  ~ . - C)

model_list <- 
  model_list %>% 
  mutate(formX = ifelse( climate_window == "none", list( formXNULL), formX  )) %>% 
  filter( !is.na(formX)) %>% 
  distinct( spp, vr, adapt_delta, climate_window, formX, left_cut )

i <- 1


for(i in 1:nrow(model_list)){ 
  
  # choose species 
  sp <- model_list$spp[i]
  ad <- model_list$adapt_delta[i]
  fx <- model_list$formX[[i]]
  lc <- model_list$left_cut[i]
  window <- model_list$climate_window[i]
  
  dat_file <- paste0('data/temp_data/', sp, '_growth_survival_dataframe.RDS')
  dat <- readRDS(dat_file)
  
  intra_comp <- paste0('W.', sp)
  
  dat$size <- scale( dat$logarea.t0 )
  
  dat$small <- as.numeric(dat$size < small)
  dat$Y    <- scale( dat$logarea.t1 )
  dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 
  dat$W.intra  <- scale( dat[ , intra_comp])
  dat$W.inter <- scale( rowSums(dat$W[, -( grep ( intra_comp , colnames(dat$W))) ] ) ) # inter specific comp. 
  
  dat <- left_censor_df(dat, left_cut = lc)
  
  if ( window != 'none' ){   
    moist <- paste0( 'C.VWC.', window)
    therm <- paste0( 'C.T.', window )
  
    # get climate effects 
    formC <- as.formula( paste0 ( '~-1 + ', moist, '*', therm  ))  ### Climate effects design matrix 
  }else if( window == 'none'){ 
    
    formC <- as.formula( '~-1')  
  }
  
  hold <- unique( dat$yid [ dat$Period == 'Modern'] ) 
  
  dl <- process_data(dat = dat, 
                       formX = fx, 
                       formC = formC,
                       formE = formE,
                       formZ = formZ, 
                       vr = vr, 
                       hold = hold, 
                     IBM = 1)
  
  mod <- rstan::stan_model(stan_model_file) # load stan model 
  
  print( paste( '### ---- species', sp, '; climate window', window, '--------------------##'))
  
  print( paste( '### ---- working on model', i, 'of', nrow(model_list),' -------------###' ))
      
  fit1 <- rstan::sampling(mod,
                            data = dl,
                            chains = nchains,
                            iter = niter,
                            cores = ncores,
                            pars = c('hold_log_lik', 'hold_SSE', 'beta', 'IBM_Y_hat'),
                            control = list(adapt_delta = ad),
                            thin = nthin, 
                            refresh = -1)
  
  saveRDS(dl, file = paste0( 'output/stan_fits/', sp, '_', vr, '_', window, '_data.RDS'))
  saveRDS(fit1, file = paste0( 'output/stan_fits/', sp, '_', vr, '_', window, '_model.RDS'))
  
}
