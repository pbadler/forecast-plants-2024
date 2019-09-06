rm(list = ls())

library(rstan)
library(tidyverse)
library(loo)

rstan_options(auto_write = T)
options( mc.cores = parallel::detectCores())

source('code/analysis/stan_data_functions.R')

testing <- T

vr <- 'growth'
model_file <- 'code/analysis/growth2_student_t.stan'

stan_model( file = model_file, 
            model_name = vr, save_dso = T)

species <- c('ARTR', 'HECO', 'POSE', 'PSSP')

if( testing ){ 
  
  k <- 2                      ### number of folds 
  n_mods <- 2
  species <- species[1:4]
  
  # STAN pars -------------- 
  ncores <- 1 
  niter <- 100
  nchains <- 1 
  nthin <- 1
  
}else{
  
  k <- 10                      ### number of folds 
  n_mods <- 8
  
  # STAN pars -------------- 
  ncores <- 4 
  niter <- 2000
  nchains <- 4 
  nthin <- 4
  
}

# --------------------------

# Model Parameters 
small <- -1                 ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(paste0 ('~ size + small + W.intra + W.inter + C')) ### Fixed effects design matrix (include climate as "C")
formE = as.formula(~ size)  ### For growth model, size dependent variance design matrix  

left_cut <- data.frame( 
  species = species, 
  lc = c(0, -1, -1, -1) 
  ) # for censored data 

# set up model selection table --------------------------# 
load('data/temp_data/climate_combos.RData')

nvars <- length(T_combos)

climate_effects <- paste0( 'C.T.', 1:nvars, '*', 'C.VWC.', 1:nvars)
climate_effects <- c('NULL', climate_effects)
climate_effects <- climate_effects[1:n_mods]

model_combos <- 
  expand.grid( species = species, model = climate_effects, fold = 1:k) %>% 
  arrange( species, model, fold) %>% 
  mutate( id = row_number() ) %>% 
  select( id, species, model, fold ) %>% 
  left_join(left_cut)

model_combos <- model_combos %>% split( .$id)
total <- length(model_combos)

dat_filenames <- dir(path = 'data/temp_data', pattern = '.*_growth_survival_dataframe.RDS', full.names = T)

prepped_dfs <- 
  mapply( x = as.character( left_cut$species), 
        y = left_cut$lc, 
        FUN = function(x, y, ... ) scale_and_fold(species = x, lc = y,...), 
        k = k, filenames = dat_filenames, SIMPLIFY = T)


counter <- 1
s <- 1
i <- 1

output <- list()

for( i in 1:length( model_combos )){ 
  
  sp <- model_combos[[i]]$species
  dat <- prepped_dfs[[sp]]  
  model <- model_combos[[i]]$model
  formC <- as.formula( paste0 ( '~-1 + ', model  ))  ### Climate effects design matrix 
  fold <- model_combos[[i]]$fold
  
  hold <- unique( dat$yid[ dat$folds == fold ] )
  
  dl <- process_data(dat = dat, 
                     formX = formX, 
                     formC = formC, 
                     formZ = formZ, 
                     formE = formE,
                     vr = vr, 
                     hold = hold, 
                     IBM = 0)
  
  print( paste( '### ---- species ', sp, ' out of ', length(species), ' -------- # '))
  print( paste( '### ---- working on rep', counter, 'of', total, ': ', 100*(counter - 1)/total, '% done ----------###' ))
  
  cat('\n\n')
  
  fit <- stan(file = model_file, 
               data = dl, 
               chains = nchains, 
               iter = niter, 
               cores = ncores,
               thin = nthin,
               pars = c('log_lik', 'hold_log_lik', 'hold_SSE'), 
               refresh = -1)
  
  stopifnot(get_num_divergent(fit) == 0) ### stop if there are divergent transitions 
  
  rhat <- range( summary( fit, 'log_lik')$summary[, 'Rhat'] )
  ins_lppd  <- sum( get_lpd(fit, 'log_lik'))
  oos_lpd  <- get_lpd(fit, 'hold_log_lik')
  oos_sse  <- summary(fit, 'hold_SSE')$summary[, 'mean']        
  hold_N   <- dl$hold_N
  N        <- dl$N 
  
  my_pars <- c('rhat', 'ins_lppd', 'oos_lpd', 'oos_sse', 'hold_N', 'N')  
  stats <- lapply( my_pars, function(x) eval(parse( text = x)))
  names( stats ) <- my_pars
  output[[i]] <- c( model_combos[[i]], stats )
   
  counter <- counter + 1 
  
  rm(list =  c(my_pars, c('fit', 'dl', 'hold', 'formC', 'model', 'dat', 'fold')))
  
}

output


