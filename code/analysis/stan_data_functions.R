extract_data <- function(df){ 
  N <- nrow(df)
  X <- df$X
  K <- ncol(X)
  Z <- df$Z 
  J <- ncol(Z)
  g <- as.numeric(factor(df$g))
  G <- length(unique(g))
  
  Y_attrib <- attributes(df$Y) # save Y center and scale for back transforming predictions

  Y <- as.numeric(df$Y)
  S <- as.numeric(df$survives)

  quad_name <- df$QuadName 
  year_name <- df$year
  
  obs <- which(df$censored == 0)
  cens <- which(df$censored == 1)
  censored <- as.numeric(df$censored)
  N_obs <- length(obs)
  N_cens <- length(cens)
  Y_obs  <- Y[obs]
  U <- unique( df$U )
  
  P1 <- df$P1
  P2 <- df$P2
  years <- unique( df$year )
  
  rm(df)
  out <- lapply( ls(), function(x) eval(parse(text = x)))
  names(out) <- ls()[!ls() == 'out']
  
  return(out)
}


split_df <- function(df, vr, hold){ 
  if(all(hold == 0) ){
    df_out <- split(df, df$g %in% hold)
    df_out$True <- data.frame(Y = rep(0,2))
    df_out$True$survives = rep( 0, 2)
    df_out$True$X = matrix(0, ncol = ncol(df$X), nrow = 2)
    df_out$True$Z = matrix(0, ncol = ncol(df$Z), nrow = 2)
    df_out$True$g = rep(0,2)
  }else if(any(hold > 0)){ 
    df_out <- split(df, df$g %in% hold)
  }
  names(df_out) <- c('train', 'hold')
  
  attr(df_out[[1]]$Y, "scaled:center") <- attr(df$Y, "scaled:center")
  attr(df_out[[1]]$Y, "scaled:scale")  <- attr(df$Y, "scaled:scale")
  
  attr(df_out[[2]]$Y, "scaled:center") <- attr(df$Y, "scaled:center") 
  attr(df_out[[2]]$Y, "scaled:scale")  <- attr(df$Y, "scaled:scale")
  
  return(df_out)
}


make_dl <- function(df){ 
  dl <- unlist( lapply( df, extract_data), recursive = F)
  names(dl) <- gsub( 'train.', '', names(dl) )
  names(dl) <- gsub( 'hold.', 'hold_', names(dl))
  return(dl)
}

get_lpd <- function(my_fit, llname = 'log_lik' ){ 
  require(loo)
  
  lpd <- log(colMeans(exp(extract_log_lik(my_fit, llname))))

  return(lpd)
}

plot_x_y <- function(myfit, X, Y, iter, bt = F){
  
  if(bt){
    Y_hat <- extract( myfit, 'Y_hat_bt')$Y_hat_bt
  }else if(!bt){
    Y_hat <- extract( myfit, 'Y_hat')$Y_hat
  }
  
  ylims = c(min(c(Y_hat[iter,], Y)), max(c(Y_hat[iter, ], Y)))
  
  par(mfrow = c(2,2))
  
  plot(X, Y, ylim = ylims)
  abline(0,1)
  
  plot(X, Y_hat[iter,], ylim = ylims)
  abline(0,1)
  
  plot(Y, Y_hat[iter, ], ylim = ylims)
  abline(0,1)
  
  par(mfrow = c(1,1))
}

left_censor_df <- function(df, left_cut){ 
  # account for left censored data 
  # left cut is the cut off for the censored data on the original scale
  # of logarea.t1
  # it is rescaled to the Y-scale and then applied 
  
  U <-  max ( df$Y[ df$logarea.t1 < left_cut ], na.rm = T )
  
  hist(df$Y)
  abline(v = U, lty = 2)
  
  df$U <- U
  df$censored <- as.numeric( df$Y <= df$U )
 
  return(df)
}


init_norm <- function(means=0, sds=1, name = NULL){
  out <- list(mapply( x = means, y = sds, function(x, y) rnorm(1, x, y)))
  names(out) <- name 
  out 
}

init_unif <- function(lower, upper, name = NULL){
  out <- list(mapply( x = lower, y = upper, function(x, y) runif(1, x, y)))
  names(out) <- name 
  out 
}


inv_logit <- function(x) { exp(x)/(1 + exp(x)) }

init_chol <- function(r_mu, r_sd, nr = 2, nc = 2, name = NULL){
  r <- 2*inv_logit(rnorm(1, r_mu, r_sd)) - 1
  out <- list(t(chol(matrix(c(1,r,r,1), 2,2))))
  names(out) <- name 
  out 
}

get_spp_and_vr <- function(dat_file, model_file){ 
  spp <- unlist( str_extract_all(c('ARTR', 'HECO', 'POSE' , 'PSSP'), string = dat_file) )
  vr <- unlist(str_extract_all(c('growth', 'survival', 'recruitment'), string = model_file))
  return(list(spp, vr))
}

process_data <- function(dat, formX, formC, formZ, formE = as.formula(~ -1), vr = 'growth', IBM = 0, ... ){
   
  if( vr == 'recruitment'){ 
    return( process_recruitment_data(dat, formX, formC, formZ, IBM = IBM, ... )) 
  }
  
  C <- model.matrix(formC, dat)
  dat$C <- scale(C)
  dat$W <- scale(dat$W)
  dat$Group <- factor(dat$gid)
  
  if( ncol(dat$C) == 0 ){ 
    formX <- update(formX, ' ~ . - C')
  }
  
  dat$X <- model.matrix(formX, data = dat)
  dat$Z <- model.matrix(formZ, data = dat)
  dat$E <- model.matrix(formE, data = dat) 
  
  dat$g <- factor(dat$yid)
  
  dat_4_IBM <- dat ### Need to preserve dataframe with NA's (dead plants) for IBM simulations 
  dat_4_IBM <- split_df(dat_4_IBM, vr = vr, hold = 0)
  dl_4_IBM <- make_dl(dat_4_IBM)
  dl_4_IBM <- dl_4_IBM[-grep('hold', names(dl_4_IBM))]
  names(dl_4_IBM) <- paste0( 'IBM_', names(dl_4_IBM))
  dl_4_IBM$IBM <- IBM
  
  if(vr == 'growth'){ 
    dat <- dat[complete.cases(dat), ]
    dat <- split_df(dat, vr, ... )
    dl <- make_dl(dat)
  }else if(vr == 'survival'){ 
    dat <- split_df(dat, vr, ... )
    dl  <- make_dl(dat)
  }
  
  return( c(dl, dl_4_IBM))
}



process_data <- function(dat, formX, formC, formZ, vr = 'growth', IBM = 0, ... ){
  
  C <- model.matrix(formC, dat)
  dat$C <- scale(C)
  dat$W <- scale(dat$W)
  dat$Group <- factor(dat$gid)
  
  if( ncol(dat$C) == 0 ){ 
    formX <- update(formX, ' ~ . - C')
  }
  
  dat$X <- model.matrix(formX, data = dat)
  dat$Z <- model.matrix(formZ, data = dat)
  
  dat$g <- factor(dat$yid)
  
  dat_4_IBM <- dat ### Need to preserve dataframe with NA's (dead plants) for IBM simulations 
  dat_4_IBM <- split_df(dat_4_IBM, vr = vr, hold = 0)
  dl_4_IBM <- make_dl(dat_4_IBM)
  dl_4_IBM <- dl_4_IBM[-grep('hold', names(dl_4_IBM))]
  names(dl_4_IBM) <- paste0( 'IBM_', names(dl_4_IBM))
  dl_4_IBM$IBM <- IBM
  
  if(vr == 'survival'){ 
    dat <- split_df(dat, vr, ... )
    dl  <- make_dl(dat)
  }else{ 
    dat <- dat[complete.cases(dat), ]
    dat <- split_df(dat, vr, ... )
    dl <- make_dl(dat)
  }
  
  return( c(dl, dl_4_IBM))
}



check_for_compiled_model <- function(vr, model_file){ 
  compiled_model <- dir('data/temp_data/', paste0(vr, '_model_compiled.RDS'), full.names = T)  
  while(length(compiled_model) == 0){
    rt <- stanc(model_file)
    stan_mod <- stan_model(stanc_ret = rt, verbose = FALSE, auto_write = F)
    saveRDS(stan_mod, file.path('data/temp_data', paste0(vr, '_model_compiled.RDS')))
    compiled_model <- dir('data/temp_data', paste0(vr, '_model_compiled.RDS'), full.names = T)
  }
  return( readRDS(compiled_model) )
}

drop_init_years <- function(my_inits, hold){ 
  
  lapply( my_inits, function(x){x$u_raw = x$u_raw[, -hold]; return(x) } )
}

drop_init_covariates <- function(my_inits, K){
  lapply( my_inits, function(x) { x$beta <- x$beta[1:K]; return(x) })  
}

find_dv_trans <- function(x){ 
  ss <-  get_sampler_params(x) 
  dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
  return(dv)
}

translate_ci <- function( ci ){ 
  lc <- (100 - ci)/2
  uc <- (100 - lc )
  return( c(lc, uc))
} 

translate_ci <- Vectorize(translate_ci)

test_pp_intervals <- function(Y_obs, my_fit, ci = c(95,90,75,50)) { 
  
  cls <- translate_ci(ci)  
  cls <- sort( as.numeric (cls/100 ))
  labs <- paste0( 'X', cls*100, '.') 
  limit_level <- factor( labs, levels = unique(labs) , labels = c(ci, rev(ci)) )
  
  labels <- data.frame( limit = labs, limit_level)
  
  temp <- summary(my_fit, 'Y_hat', c(cls, 0.5) )$summary
  
  data.frame( temp, Y = Y_obs ) %>%
    mutate( id = row_number()) %>% 
    rename( 'med' = `X50.` ) %>% 
    gather( limit, value, starts_with('X')) %>% 
    left_join(labels) %>% 
    mutate( side = ifelse( value < med, 'lower', 'upper' )) %>% 
    select( - limit ) %>% 
    spread( side, value) %>%  
    mutate( below = Y < lower, 
            above = Y > upper )
} 

scale_and_fold <- function(species, vr, lc = 0, small = -1,  k = 10, filenames){
  
  file_name <- filenames[ str_detect(filenames, species)  ]
  
  dat <- readRDS(file_name)
  dat <- dat[ dat$Period == 'Historical',  ] 
  
  folds <- kfold_split_grouped(k, dat$yid) 
  dat$folds <- folds
  
  k_folds <- 
    dat %>% 
    distinct(yid, folds)
  
  if(vr != 'recruitment'){
    dat$size <- scale( dat$logarea.t0 )
    dat$small <- as.numeric(dat$size < small)
    dat$Y    <- scale( dat$logarea.t1 )
    intra_comp <- paste0('W.', species)
    dat$W.intra  <- scale( dat[ , intra_comp])
    dat$W.inter <- scale( rowSums(dat$W[, -( grep ( intra_comp , colnames(dat$W))) ] ) ) # inter specific comp. 
  
    dat <- left_censor_df(dat, lc )  
    
  }else if( vr == 'recruitment'){ 
    
    dat <- 
      dat %>% 
      ungroup %>% 
      rowwise( ) %>% 
      mutate( total_basal_cover = cov.HECO + cov.POSE + cov.PSSP) %>% 
      mutate( total_open = 100*100 - total_basal_cover) %>% 
      mutate( l_open = log(total_open))
    
    dat$W <- scale( dat$total_basal_cover)
    
    dat <- as.data.frame( dat ) 
  }
  
  dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 
  
  return(dat)
}


get_model_score <- function( model_combos,
                             vr,
                             model_file, 
                             prepped_dfs, 
                             formZ = as.formula(~ size ), 
                             formX = as.formula(~ size + small + W.intra + W.inter + C),
                             n_iter = 2000, 
                             nthin = 4){ 
  id <- model_combos$id                      
  sp <- model_combos$species
  ad <- model_combos$ad
  model <- model_combos$model
  formC <- as.formula( paste0 ( '~-1 + ', model  ))  ### Climate effects design matrix 
  fold <- model_combos$fold
  
  dat <- prepped_dfs[[sp]]  
  hold <- unique( dat$yid[ dat$folds == fold ] )
  

  dl <- process_data(dat = dat, 
                     formX = formX, 
                     formC = formC, 
                     formZ = formZ, 
                     formE = ~ 1,
                     vr = vr, 
                     hold = hold, 
                     IBM = 0)
  
  if( vr == 'growth'){ 
    ad <- 0.8
  }
  
  fit <- stan(file = model_file, 
              data = dl, 
              chains = 4, 
              iter = n_iter, 
              cores = 4,
              thin = nthin,
              control = list(adapt_delta = ad), 
              verbose = F,
              pars = c('log_lik', 'hold_log_lik', 'hold_SSE'), 
              refresh = 0)
  
  nd <- get_num_divergent( fit)
  if( nd > 0 ) stop( paste0( nd,  ' Divergent transitions in model run ', id ))
  
  if(dim( summary(fit)$c_summary)[3] == 4){ 
    stop( paste('missing chains for run', id ) )
  }
  
  mx_rhat <- max(summary( fit, 'log_lik')$summary[, 'Rhat'])
  mn_rhat <- min( summary( fit, 'log_lik')$summary[, 'Rhat'])
  ins_lppd  <- sum( get_lpd(fit, 'log_lik'))
  oos_lpd  <- sum( get_lpd(fit, 'hold_log_lik'))
  oos_sse  <- summary(fit, 'hold_SSE')$summary[, 'mean']        
  hold_N   <- dl$hold_N
  N        <- dl$N 
  
  my_pars <- c('mx_rhat', 'mn_rhat', 'ins_lppd', 'oos_lpd', 'oos_sse', 'hold_N', 'N')  
  stats <- lapply( my_pars, function(x) eval(parse( text = x)))
  names( stats ) <- my_pars
  
  temp_scores <- c( model_combos, stats, vr = vr )
  outfile <- paste0('output/model_scores/', vr, '_scores_', id, '.rds')
  write_rds( temp_scores, path = outfile) 
}

