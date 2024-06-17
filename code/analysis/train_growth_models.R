rm(list = ls() )
library(tidyverse)
library(lubridate)
library(glmnet)
source('code/analysis/functions.R')
# For each species this script fits a: 
#  1. Growth model with size x climate interaction. 
#  2. Null growth model without climate 
#  3. Growth model for small plants 
#  4. Null growth model for small plants 
#  
#  Finally it saves training and testing data for each model.
# ----------------------------------------------------------  

split_year <- 2010 # Training testing split year 
size_cutoff <- -1  # log scale size cutoff between large and small 
quad_info <- read_csv( file = 'data/quad_info.csv')
daily_weather <- read_csv('data/temp/daily_weather_for_models.csv')

### make seasonal_weather dataframe  

daily_weather$year <- year(daily_weather$date)

seasonal_weather <- daily_weather %>%
  group_by(year,season,Treatment) %>%
  summarise(VWC = mean(VWC), TAVG = mean(TAVG)) 

# associate spring with year t, not t + 1
seasonal_weather$year[seasonal_weather$season=="spring"] <- seasonal_weather$year[seasonal_weather$season=="spring"] - 1

# change to wide format
seasonal_weather <- pivot_wider(seasonal_weather, names_from=season, values_from = c(VWC,TAVG))

# add lagged weather
tmp <- seasonal_weather
tmp$year <- tmp$year + 1
names(tmp)[3:ncol(tmp)] <- paste0(names(tmp)[3:ncol(tmp)],"_lag")
seasonal_weather <- merge(seasonal_weather,tmp)
rm(tmp)

### loop through species

species_list <- c('HECO')# ,'HECO', 'POSE', 'PSSP')

for( sp in species_list){ 
  
  W.intra <- paste0( 'W.', sp)
  
  # import demography data
  size <- read.csv(paste0("data/temp/",sp,"_growth.csv"))
  
  # join to climate data
  size <- merge(size, seasonal_weather,all.x=T)
  
  # split by size
  size_small <- subset(size, logarea.t0 <= size_cutoff )
  size <- subset(size, logarea.t0 > size_cutoff )
  
  # prepare data
  
  # make model matrix
  X <- data.frame(size = size$logarea.t0,
                  intra = size[,which(names(size)==W.intra)],
                  VWC_fall = size$VWC_fall,
                  VWC_spring = size$VWC_spring,
                  VWC_spring_lag = size$VWC_spring_lag,
                  T_fall = size$TAVG_fall,
                  T_winter = size$TAVG_winter,
                  T_spring = size$TAVG_spring,
                  T_spring_lag = size$TAVG_spring_lag,
                  sizeVWC_fall = size$logarea.t0*size$VWC_fall,
                  sizeVWC_spring = size$logarea.t0*size$VWC_spring,
                  sizeVWC_spring_lag = size$logarea.t0*size$VWC_spring_lag,
                  sizeT_fall = size$logarea.t0*size$TAVG_fall,
                  sizeT_winter = size$logarea.t0*size$TAVG_fall,
                  sizeT_spring = size$logarea.t0*size$TAVG_spring,
                  sizeT_spring_lag = size$logarea.t0*size$TAVG_spring_lag)
  X <- as.data.frame(scale(X))  # standardize the continuous covariates
 
  # split training and testing 
  X_train <- X[size$year <= split_year,]
  X_test <-  X[size$year > split_year,]
  y_train <- size$logarea.t1[size$year <= split_year]
  y_test <- size$logarea.t1[size$year > split_year]
  rm(X)
  
  ### fit models
  
  # Null model  
  mod_null <- lm( y_train ~ size + intra, data=X_train)
  model_name <- paste0( sp, '_growth_null')
  saveRDS(mod_null, paste0( "output/growth_models/", model_name, ".rds"))
  
  # Climate informed model 
  
  # don't penalize size or previous year's abundance
  pen_facts <- c(0,0,rep(1,ncol(X_train)-2))
  
  # specify leave-one-year-out cross-validation
  my_folds <- as.numeric(as.factor(size$year[size$year <= split_year]))
  
  # Fit the model at all levels of lambda
  fit <- cv.glmnet(
    x = as.matrix(X_train),
    y = y_train, 
    family = "gaussian",
    alpha = 0,  # 0 == ridge regression, 1 == lasso, 0.5 ~~ elastic net
    type.measure="mae",
    penalty.factor = pen_facts,
    foldid = my_folds
  )
  
  # look at CV score vs penalty plot
  plot(log(fit$lambda),fit$cvm)
  
  # extract and look at the best coefficients
  best_coefs = fit$glmnet.fit$beta[,which(fit$lambda==fit$lambda.min)]
  print(best_coefs)
  
  # plot coefficients against lambda
  matplot(log(fit$lambda),t(fit$glmnet.fit$beta),ylab="Estimates",type="l")
  abline(v=log(fit$lambda.min),col="red")
  
  # compare predictions for training data
  null_pred <- predict(mod_null) 
  null_mae <- mean(abs(y_train - null_pred))
  clim_pred <- predict(fit, newx=as.matrix(X_train), s="lambda.min")
  clim_mae <- mean(abs(y_train - clim_pred))
  print(c(null_mae,clim_mae))
  
  ### make predictions for test data
  null_pred <- predict(mod_null,new=X_test) 
  null_mae <- mean(abs(y_test - null_pred))
  clim_pred <- predict(fit,newx = as.matrix(X_test),s="lambda.min")
  clim_mae <- mean(abs(y_test - clim_pred))
  print(c(null_mae,clim_mae))
  
}
  
  
  
  
  
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   saveRDS(my_mod, paste0( "output/growth_models/", model_name, ".rds"))
# 
#   # Null model  
#   frm_null <- as.formula(  frm <- as.formula( paste0( 'area ~ area0 + ', W.intra,   ' + (1|year/Group)' )))
#   my_mod_null <- lmer( formula = frm_null, data = training_large, control = control_lmer, REML = T)
#   model_name <- paste0( sp, '_growth_null')
#   saveRDS(my_mod_null, paste0( "output/growth_models/", model_name, ".rds"))
#   
#   # no intxn models 
#   temp <- 
#     growth_windows %>% 
#     filter( species == sp, 
#             type == 'growth_no_intxn')
#   
#   first_var <- temp[ temp$fit == 1, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')] 
#   second_var <- temp[ temp$fit == 2, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')]
#   
#   clim1 <- getWindowAvg(daily_weather, 
#                         var = str_remove( first_var$climate[1], '_scaled'), 
#                         open = first_var$WindowOpen, 
#                         close = first_var$WindowClose)
#   
#   clim2 <- getWindowAvg(daily_weather, 
#                         var = str_remove( second_var$climate[1], '_scaled'), 
#                         open = second_var$WindowOpen, 
#                         close = second_var$WindowClose)
#   
#   
#   temp_clim <- clim1 %>% 
#     left_join(clim2, by = c('Treatment', 'year')) %>% 
#     mutate( !!first_var$climate[1]:= as.numeric( scale(.data[[str_remove(first_var$climate[1], '_scaled')]] ))) %>% 
#     mutate( !!second_var$climate[1]:=as.numeric( scale(.data[[str_remove(second_var$climate[1], '_scaled')]]))) %>% 
#     dplyr::select( Treatment, year, 
#                  all_of(str_remove( first_var$climate[1], '_scaled')), 
#                  all_of(str_remove( second_var$climate[1], '_scaled')), 
#                  all_of(first_var$climate[1]), 
#                  all_of(second_var$climate[1]))
#   
#   
#   # split training and testing ------------- # 
#   temp_dat <- 
#     size %>% 
#     filter( size_class == 'large') %>% 
#     filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>% 
#     ungroup() %>%
#     left_join(temp_clim, by = c('year', 'Treatment')) %>% 
#     mutate( Split = ifelse( year <= split_year, 'Training', 'Testing')) %>% 
#     split( .$Split)  
#   
#   temp_name <- paste0( sp, '_growth_no_intxn')
#   
#   write_csv( temp_dat$Training, paste0( 'data/temp/', temp_name, '_training_data.csv'))
#   write_csv( temp_dat$Testing, paste0( 'data/temp/', temp_name, '_testing_data.csv'))
#   
#   # --------------------------------------------------------- #
#   training_no_intxn <- temp_dat$Training %>% filter( size_class == 'large')
#   
#   # Climate informed model 
#   frm <- as.formula( paste0( 'area ~ ', 'area0  +' , first_var$climate[1], ' + ', second_var$climate[1], " + " , W.intra,   ' + (1|year/Group)' ))
#   
#   my_mod <- lmer( formula = frm, data = training_no_intxn, control = control_lmer, REML = T)
#   saveRDS(my_mod, paste0( "output/growth_models/", temp_name, ".rds"))
#   
#   # small plant models 
#   temp <- 
#     growth_windows %>% 
#     filter( species == sp, 
#             type == 'growth_small')
#   
#   first_var <- temp[ temp$fit == 1, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')] 
#   second_var <- temp[ temp$fit == 2, c('climate', 'WindowOpen', 'WindowClose', 'index', 'fit')]
#   
#   clim1 <- getWindowAvg(daily_weather, 
#                         var = str_remove( first_var$climate[1], '_scaled'), 
#                         open = first_var$WindowOpen, 
#                         close = first_var$WindowClose)
#   
#   clim2 <- getWindowAvg(daily_weather, 
#                         var = str_remove( second_var$climate[1], '_scaled'), 
#                         open = second_var$WindowOpen, 
#                         close = second_var$WindowClose)
#   
#   temp_clim <- clim1 %>% 
#     left_join(clim2, by = c('Treatment', 'year')) %>% 
#     mutate( !!first_var$climate[1]:= as.numeric( scale(.data[[str_remove(first_var$climate[1], '_scaled')]] ))) %>% 
#     mutate( !!second_var$climate[1]:=as.numeric( scale(.data[[str_remove(second_var$climate[1], '_scaled')]]))) %>% 
#     dplyr::select( Treatment, year, 
#                  all_of(str_remove( first_var$climate[1], '_scaled')), 
#                  all_of(str_remove( second_var$climate[1], '_scaled')), 
#                  all_of(first_var$climate[1]), 
#                  all_of(second_var$climate[1]))
#   
#   
#   # split training and testing ------------- # 
#   temp_dat <- 
#     size %>% 
#     filter( size_class == 'small') %>% 
#     filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>% 
#     ungroup() %>%
#     left_join(temp_clim, by = c('year', 'Treatment')) %>% 
#     mutate( Split = ifelse( year <= split_year, 'Training', 'Testing')) %>% 
#     split( .$Split)  
#   
#   temp_name <- paste0( sp, '_growth_small')
#   
#   write_csv( temp_dat$Training, paste0( 'data/temp/', temp_name, '_training_data.csv'))
#   write_csv( temp_dat$Testing, paste0( 'data/temp/', temp_name, '_testing_data.csv'))
#   
#   # --------------------------------------------------------- #
#   training_small <- temp_dat$Training %>% filter( size_class == 'small')
#   
#   # Climate informed model 
#   frm <- as.formula( paste0( 'area ~ ', first_var$climate[1], ' + ', second_var$climate[1], " + " , W.intra,   ' + (1|year)' ))
#   my_mod <- lmer( formula = frm, data = training_small, control = control_lmer, REML = T)
#   saveRDS(my_mod, paste0( "output/growth_models/", temp_name, ".rds"))
#   
#   # Null model  
#   frm_null <- as.formula(  frm <- as.formula( paste0( 'area ~ ', W.intra,   ' + (1|year)' )))
#   my_mod_null <- lmer( formula = frm_null, data = training_small, control = control_lmer, REML = T)
#   
#   saveRDS(my_mod_null, paste0( "output/growth_models/", temp_name, "_null.rds"))
#   
#     
# }
