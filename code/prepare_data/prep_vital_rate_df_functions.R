# clean dataframe functions 
# functions for preparing vital rate and climate data for analysis 
# scale climate data
# merge growth and survival 
# merge vital rate with climate data

# functions -------------------------------------------------------------------------- # 

set_up_basic_vars <- function(x){ 
  
  require(stringr)
  
  x$Treatment <- factor(x$Treatment)
  x$obs_id <- as.numeric( row.names(x))
  x$yid <- as.numeric( factor(x$year)) 
  x$gid <- as.numeric( factor(x$Group))
  x$quad <- as.numeric(str_extract(x$quad, '[0-9]+$')) 
  x$treat <- as.numeric(factor(x$Treatment))
  x$spp <- as.numeric(factor( x$species))
  
  return(x)
}

clean_growth_survival <- function(clim, growth_file, survival_file){
  gdat <- readRDS( growth_file)
  sdat <- readRDS( survival_file)
  
  # merge growth, survival and climate data -------------------------------------- #
  sdat$logarea.t0 <- sdat$logarea
  x <- merge(sdat, gdat, all.x = T)
  x <- merge(x, clim)
  
  x$Period <- ifelse(x$year > 2010, "Modern", "Historical")
  
  x <- set_up_basic_vars(x)
  
  #x$X <- scale(x$logarea.t0) # center X on mean.  Center BEFORE splitting up survival and growth, modern and historical
  
  W <- x[ , grep ( '^W\\.', names( x))]
  W <- as.matrix( W )[,1:6] # big four competition effects
  W <- scale( W ) # scale competition
  x$W <- W
  return(x)
}


clean_recruitment <- function(clim, recruitment_file){ 
  
  require(stringr)
  
  x <- readRDS( recruitment_file)
  x$Period <- ifelse(x$year > 2010, "Modern", "Historical")
  x <- merge(x, clim )
  
  x <- set_up_basic_vars(x) 
  
  #x$gm <- model.matrix.lm(~ x$Group)
  #x$tm <- model.matrix.lm(~ x$Treatment)[ , -1]
  
  #x$C <- as.matrix( x[ ,  grep( '^C\\.', names(x) ) ] )
  #colnames(x$C) <- str_replace( colnames( x$C ) , '^C\\.' , '')
  
  x$parents1 <- as.matrix( x[ , grep( '^cov\\.', names(x))] )
  x$parents2 <- as.matrix( x[ , grep( '^Gcov\\.', names(x))] )
  
  return(x)
}



make_my_climate_combos <- function(prefix = 'C.VWC.'){
  vars1 <- paste0(prefix, c('sp.l', 'sp.0', 'sp.1'))
  vars2 <- paste0(prefix, c('f.1', 'f.0', 'su.1', 'su.0'))
  
  c1 <- lapply( 1:3, combn, x = vars1, simplify = F)
  c1 <- unlist( c1 , recursive = F)
  c1 <- c1[unlist(lapply( c1, function(x, reqs) all( reqs %in% x ), reqs = c(vars1[3])))]
  c1 <- c1[-2]
  
  c2 <- lapply(1:4, combn, x = vars2, simplify = F)
  c2 <- unlist(c2, recursive = F)
  c2 <- c2[ grep('su.1', c2) ] 
  c2 <- c2[ -c(3, 5, 6, 7)]
  
  
  list(
    c(c1[[1]]), 
    c(c1[[2]]), 
    c(c1[[3]]),
    c(c1[[1]], c2[[1]]),
    c(c1[[1]], c2[[2]]), 
    c(c1[[3]], c2[[3]]), 
    c(c1[[3]], c2[[4]]))
  
}

add_climate_combos <- function(dat){
  
  VWC_combos <- make_my_climate_combos()
  T_combos <- make_my_climate_combos(prefix = 'C.T.')
  
  save(VWC_combos, T_combos, file = 'data/temp_data/climate_combos.RData')
  
  VWC.spring <- lapply( VWC_combos, function(x,df) apply( as.matrix( df[,x] ), 1, mean ) , df = dat)
  T.spring <- lapply( T_combos, function(x,df) apply( as.matrix(df[,x]), 1, mean), df = dat)
  
  names(VWC.spring) <- paste0('C.VWC.', 1:length(VWC.spring))
  names(T.spring) <- paste0('C.T.', 1:length(T.spring))
  
  cbind(dat, cbind( do.call(cbind, VWC.spring), do.call(cbind, T.spring)))
}

scale_climate_vars <- function( clim_file, clim_vars){ 
  clim <- readRDS(clim_file)  
  clim <- clim[ ,c('Treatment', 'year', clim_vars)]
  names ( clim ) [ grep( '^(P\\.)|(VWC\\.)|(T\\.)',  names( clim ) ) ] <- paste0 ( 'C.', names( clim ) [ grep( '^(P\\.)|(VWC\\.)|(T\\.)', names(clim ) ) ] )
  clim <- clim[ complete.cases(clim), ] 
  clim <- add_climate_combos( clim )
  clim[ , grep( 'C\\.', names(clim))] <- scale(clim[ , grep( 'C\\.', names(clim))]) # scale climate data 
  
  return(clim)
}

