library(lme4)
library(optimx)
library(dfoptim)
# LMER optimization options
control_lmer = lmerControl(
  optimizer = "optimx",
  calc.derivs = FALSE,
  optCtrl = list(
    method = "nlminb",
    starttests = FALSE,
    kkt = FALSE
  )
)
control_lmer$optCtrl$eval.max <- 1e8
control_lmer$optCtrl$iter.max <- 1e8


prep_growth_for_climWin <- function( species, last_year, quad_info, size_cutoff = -1 ) { 
  
  intra <- paste0( 'W.' , species )
  dat <- read_csv( paste0 ( 'data/temp/', paste0( species, '_size.csv')))
  dat <- dat[ dat$year < last_year, ]
  mnyear <- min(dat$year)
  mxyear <- max(dat$year)
  pid <- unique( dat$pid )
  
  dat <- 
    expand.grid(year = mnyear:mxyear, pid = unique(pid)) %>% 
    mutate( quad = str_extract(pid, "[A-z0-9]+(?=_)"))   %>% 
    left_join(dat, by = c('quad', 'pid', 'year'))  %>% 
    arrange(pid, year) 
  
  my_dat <- 
    dat %>% 
    dplyr::select(quad, pid, age, year, area, starts_with('W.'))  %>% 
    mutate( W.intra = eval(parse(text = intra)))
  
  growth <- my_dat %>% 
    left_join(quad_info, by = 'quad') %>% 
    group_by( pid ) %>% 
    arrange(pid, year) %>% 
    mutate(area = log(area))  %>% 
    mutate( area0 = lag(area)) %>% 
    arrange( pid, year ) %>% 
    filter(!is.na( area ), !is.na(area0))  %>% 
    filter( area0 > size_cutoff ) %>% 
    mutate( climate = 1) %>% 
    mutate( date_reformat  = paste0( "15/06/", year)) 
  
  return( growth )
}

prep_survival_for_climWin <- function( species, last_year, quad_info, size_cutoff = -1 ) { 
  
  intra <- paste0( 'W.' , species )
  dat <- read_csv( paste0 ( 'data/temp/', paste0( species, '_survival.csv')))
  dat <- dat[ dat$year < last_year, ]
  mnyear <- min(dat$year)
  mxyear <- max(dat$year)
  pid <- unique( dat$pid )
  
  dat <- 
    expand.grid(year = mnyear:mxyear, pid = unique(pid)) %>% 
    mutate( quad = str_extract(pid, "[A-z0-9]+(?=_)"))   %>% 
    left_join(dat, by = c('quad', 'pid', 'year'))  %>% 
    arrange(pid, year) 
  
  my_dat <- 
    dat %>% 
    dplyr::select(quad, pid, age, year, area, starts_with('W.'))  %>% 
    mutate( W.intra = eval(parse(text = intra)))
  
  growth <- my_dat %>% 
    left_join(quad_info, by = 'quad') %>% 
    group_by( pid ) %>% 
    arrange(pid, year) %>% 
    mutate(area = log(area))  %>% 
    mutate( area0 = lag(area)) %>% 
    arrange( pid, year ) %>% 
    filter(!is.na( area ), !is.na(area0))  %>% 
    filter( area0 > size_cutoff ) %>% 
    mutate( climate = 1) %>% 
    mutate( date_reformat  = paste0( "15/06/", year)) 
  
  return( growth )
}


prep_survival_for_climWin <- function( species, last_year, quad_info){ 
  
  intra <- paste0( 'W.' , species )
  
  dat <- read_csv( paste0 ( 'data/temp/', paste0( species, '_survival.csv'))) %>% 
    mutate( pid = paste( quad, trackID, sep = "_"))
  
  dat <- dat %>% 
    mutate( year = year + 1) %>% # Add one year so that the survival "response" coincides with the current year, not the future
    filter( year < last_year )  
  
  mnyear <- min(dat$year)
  mxyear <- max(dat$year)
  pid <- unique( dat$pid )
  
  my_dat <- 
    dat %>% 
    dplyr::select(quad, pid, age, year, logarea, survives, starts_with('W.'))  %>% 
    mutate(W.intra = eval(parse(text = intra))) %>% 
    rename( "area0" = logarea) # for consistency with growth data 
  
  survival <- 
    my_dat %>% 
    group_by( pid ) %>% 
    arrange(pid, year) %>% 
    filter( !is.na(area0), 
            !is.na(survives),
            is.finite(area0), 
            year > 1925) %>%
    mutate( date_reformat = paste0( '15/06/', year)) %>% 
    left_join(quad_info, by = 'quad') %>% 
    mutate( climate = 1 )
  
  return( survival )
}

addVars <- function( climWin , data1 , responseVar = 'survives'){ 
  
  bestModel <- climWin$combos %>% filter( DeltaAICc == min(DeltaAICc))
  bestVar <- bestModel$climate %>% as.character()
  bestModelData <- climWin[[ which( climWin$combos$climate == bestVar )]]
  
  data2 <- bestModelData$BestModelData %>% 
    rename(  !!bestVar:=climate ) %>% 
    rename(  !!responseVar:= yvar )  %>% 
    mutate( climate = 1 )
  
  data2$date_reformat <- data1$date_reformat
  
  if( bestVar == 'VWC_scaled'){ 
    addVars = c('TMAX_scaled', 'TAVG_scaled', 'TMIN_scaled')
  }else{ 
    addVars = 'VWC_scaled'  
  }

  return( list(data2 = data2, bestVar = bestVar, addVars = addVars ))
}

# Functions for processing slidingwin results: 
get_calendar_dates <- function( x , foo_year = 1999, ref_date = '06-15') { 
  # x is results combos from sliding win 
  foo_ref_date <- paste( foo_year, ref_date, sep = '-')
  
  x %>% 
    mutate( foo_date_open = ymd(foo_ref_date) - WindowOpen*7,  
            foo_date_close = ymd(foo_ref_date) - WindowClose*7) %>% 
    mutate( month_day_open = paste( month.abb[month(foo_date_open)], day(foo_date_open), sep = '-'), 
            month_day_close = paste( month.abb[month(foo_date_close)], day(foo_date_close), sep = '-'), 
            years_back_open = floor( foo_year - year( foo_date_open )), 
            years_back_close = floor( foo_year - year(foo_date_close))) %>% 
    mutate( Open = paste( month_day_open, years_back_open, 'year back'), 
            Close = paste( month_day_close, years_back_close, 'year back')) %>% 
    mutate( Open = str_replace(Open, '0 year back', 'current year'), 
            Close = str_replace(Close, '0 year back', 'current year')) %>% 
    select( -starts_with('foo'), -starts_with('month'), -starts_with('year'))
}

my_delta_plot <- function( climWinFit, temp_par ){ 
  
  delta_plot <- plotdelta(climWinFit[[temp_par$fit]][[temp_par$index]]$Dataset)
  
  delta_plot <- 
    delta_plot + 
    guides(fill = guide_colorbar( 'DeltaAICc')) + 
    ggtitle( 
      paste0(
        LETTERS[temp_par$fit] , ') ', 
        paste( unlist( temp_par[, 'climate']), collapse = ', ')
        )) + 
    theme( plot.title = element_text( size = 10, hjust = 0), 
           legend.position = 'right', 
           axis.text =  element_text(face = 'plain'), 
           axis.title =  element_text(face = 'plain'))
}

size_by_climate_plots <- function( model, clim_var, label = 'D') { 
  
  frame <- model@frame 
  type <- 'survival'
  if(class(model) == 'lmerMod'){ 
    model <- update(model, data = model@frame, REML = T)
    type <- 'growth'
  }
  W.intra <- median( frame$W.intra, na.rm = T)
  
  mgrid <- expand.grid( climate = seq( min(frame$climate), max(frame$climate), length.out = 100), 
                        area0 = seq( min(frame$area0), max(frame$area0), length.out = 100))
  
  v2name <- names(frame)[ ! names(frame)  %in% names(mgrid)]
  v2name <- v2name[ str_detect( v2name, pattern = 'TMAX|TMIN|TAVG|VWC') ]
  
  title <- paste0( label, ') ', ' size by ', str_extract(clim_var, '[A-Z]+'))
  
  if(length(v2name) == 1 ){ 
    v2 = median( frame %>% pull( !!v2name), na.rm = T)
    
    mgrid <- mgrid %>% 
      mutate( !!v2name:= v2 )
    
    title <- paste0( label, ') ', ' size by ', str_extract(clim_var, '[A-Z]+'), ' with ', v2name, '@', round( v2, 2) )
  } 
  
  mgrid$W.intra <- W.intra 
  mgrid$yhat <- predict( model, newdata = mgrid , re.form = NA)
  
  
  if( type == 'growth'){ 
    
    mgrid %>% 
      mutate( growth = yhat - area0 ) %>% 
      ggplot( aes( x = area0, y = climate, z = growth ))  + 
      geom_tile( aes( fill = growth )) + 
      geom_contour( ) + 
      ylab(clim_var ) + 
      scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') + 
      ggtitle(title ) + 
      theme( plot.title = element_text( size = 8, hjust = 0))
    
  }else if( type == 'survival'){ 
    
    mgrid %>% 
      rename( survival = yhat) %>% 
      ggplot( aes( x = area0, y = climate, z = survival ))  + 
      geom_tile( aes( fill = survival )) + 
      geom_contour( ) + 
      ylab(clim_var ) + 
      scale_fill_gradient2('log odds ratio', low = 'blue', mid = 'white', high = 'red') + 
      ggtitle(title ) + 
      theme( plot.title = element_text( size = 8, hjust = 0))
  }
  
}

getBestWindows <- function( combos ){ 
  topModel <- which.min(combos$DeltaAICc)
  window <- combos[topModel, ] %>% 
    mutate(index = topModel)
  return( window)
}

makeWindowTable <- function( WindowResult, species_name, type_name ) { 
  
  getBestWindows(WindowResult[[1]]$combos) %>% 
    mutate( fit = 1 ) %>% 
    bind_rows(
      getBestWindows(WindowResult[[2]]$combos) %>% 
        mutate( fit = 2 )
    ) %>% 
    mutate( species = species_name, type = type_name) %>% 
    mutate( climate = as.character( climate))
}


make_window_plots <- function( window_df, sp ){ 
  
  temp_windows <- 
    window_df %>% 
    filter( species == sp)
  
  temp_par <- temp_windows[1, c('species', 'type', 'climate', 'fit', 'index', 'WindowOpen', 'WindowClose')]
  
  # load results file for the current species
  if( temp_par$type[1] == 'growth'){ 
    load(file = paste0( 'output/growth_models/', sp,  '_', temp_par$type[1], '_mer_monthly_ClimWin.rda'))
    
  }else if( temp_par$type[1] == 'survival'){ 
    load(file = paste0( 'output/survival_models/', sp,  '_', temp_par$type[1], '_mer_monthly_ClimWin.rda'))
  }
  temp_obj <- ls()[str_detect( ls(), sp) & str_detect(ls(), temp_par$type[1])]  
  results <- eval(parse(text = temp_obj))  
  
  # First climate variable ----------------# 
  fit1  <- my_delta_plot(results, temp_par)
  
  sXc1 <- size_by_climate_plots(
    results[[temp_par$fit]][[temp_par$index]]$BestModel, 
    clim_var = paste0( temp_par$climate, ' (', temp_par$WindowOpen, '-', temp_par$WindowClose, ' mos. back)')) 
  
  # Second fit, second climate variable ---# 
  temp_par <- temp_windows[2, c('species', 'type', 'climate', 'fit', 'index', 'WindowOpen', 'WindowClose')]
  fit2  <- my_delta_plot(results, temp_par)
  
  sXc2 <- size_by_climate_plots(
    results[[temp_par$fit]][[temp_par$index]]$BestModel, 
    clim_var = paste0( temp_par$climate, ' (', temp_par$WindowOpen, '-', temp_par$WindowClose, ' mos. back)'), 
    label = 'E') 
  
  temp_title <- paste0( temp_par[, c('species', 'type')], collapse = ' ')
  
  png(filename = paste0('figures/', str_squish(temp_title), "_window_plots.png"),
      width = 7,
      height = 6,
      units = 'in',
      res = 300)

  grid.arrange(fit1, fit2, sXc1, sXc2,
               top = text_grob(temp_title, size = 12))
  dev.off() 
}


getWindowAvg <- function(weather, var, open, close, ref_day = '06-15'){ 
  # Intended to recreate ClimWin Windows 
  # very close for most cases but slightly different avg.
  # with windows more than two months.  
  
  weather$ref_day <- mdy( paste( ref_day, weather$year , sep = '-'))
  
  windows <- weather %>% 
    distinct(Treatment, ref_day) %>% 
    mutate( open_date = rollback( date(ref_day - duration( open, 'months') ), preserve_hms = T)) %>% 
    mutate( close_date = rollforward( date(ref_day - duration( close, 'months')), preserve_hms = T)) %>% 
    mutate( year = year( ref_day))
  
  climate <- NA
  x <- unlist( weather[, var] )
  
  for( i in 1:nrow( windows)){ 
    climate[i] <- mean( x[ weather$Treatment == windows$Treatment[i] & weather$date > windows$open_date[i] & weather$date <= windows$close_date[i] ])
  }
  
  windows %>% 
    mutate( !!var:= climate )
  
}

cross_validate_growth <- function ( model, data, folds){ 
  out <- list( R2 = NA, RMSE = NA, MAE = NA)
  
  for ( f in folds ) { 
    train = data[ -f, ]  
    test = data[ f, ]
    
    temp_model <- update( model,  data = train, REML = T) 
    yhat <- predict( model, newdata = test , allow.new.levels = T)
    
    out$R2 <- c(out$R2, caret::R2( pred = yhat, obs = test$area ))
    out$RMSE <-c(out$RMSE,  caret::RMSE(pred = yhat, obs = test$area ))
    out$MAE <- c(out$MAE, caret::MAE( pred = yhat, obs = test$area))
  }
  
  return( do.call( cbind, out )[-1, ] %>% data.frame() )
}


cross_validate_survival <- function ( model, data, folds ){ 
  out <- list( AUC = NA , RMSE = NA )
  
  for ( f in folds ) { 
    train = data[ -f, ]  
    test = data[ f, ]

    temp_model <- update( model,  data = train) 
    pred <- predict( model, newdata = test , allow.new.levels = T, type = 'response')
    
    if ( n_distinct(test$survives) == 2 ){ 
      tempROC <-  pROC::roc(test$survives, pred)
      tempAUC <- pROC::auc(tempROC) 
    }else if( n_distinct(test$survives) != 2){ 
      print( "single level survival data")
      tempROC <- NA
      tempAUC <- NA
    }
    
    
    out$AUC <- c(out$AUC, tempAUC)
    out$RMSE <-c(out$RMSE,  caret::RMSE(pred = pred, obs = test$survives ))
    
  }
  
  return(  do.call( cbind, out)[-1 , ] %>% data.frame() )
  
}




