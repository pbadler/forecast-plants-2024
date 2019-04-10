rm(list = ls())
library(tidyverse)

stan_mods <- read_csv('~/Dropbox/projects/ExperimentTests/precip/output/model_ranks_new.csv')

top_mods <- 
  stan_mods %>% 
  group_by( vr, spp) %>% 
  filter( oos_lppd == max(oos_lppd) ) %>% 
  select( vr, spp, climate_window) 

out <- list(NA)

for( i in 1:nrow(top_mods)){ 
  vr <- top_mods$vr[i]
  spp <- top_mods$spp[i]
  window <- top_mods$climate_window[i]
  
  
  if( !is.na(window) ){
    fn <- paste0( 'output/stan_fits/', spp, '_', vr, '_model_', window, '_top_model.RDS')
    fit <- readRDS(fn)
    beta <- rstan::summary(fit, c('beta[6]', 'beta[7]', 'beta[8]'))$summary[, c('mean', '2.5%', '50%', '97.5%')]
    beta <- 
      data.frame( t(beta) ) %>% 
      mutate( stat = row.names(.))
    
    out[[i]] <- data.frame( spp  = spp, vr = vr, climate_window = window, beta ) 
  }
}

top_mods <- 
  top_mods %>% 
  left_join( do.call(rbind, out ) ) %>% 
  rename( 'Moist' = beta.6., 'Temp' = beta.7., 'Temp_x_Moist' = beta.8.) %>% 
  gather( par, est, Moist:Temp_x_Moist) %>% 
  spread(stat, est) %>%
  select( - `<NA>`) 

gg1 <- 
  top_mods %>% 
  ggplot(aes( x = par, y = `50%`, ymin = `2.5%`, ymax = `97.5%`)) + 
  geom_point() + 
  geom_errorbar() + 
  geom_hline(aes( yintercept = 0), linetype = 2) + 
  facet_grid(vr~spp, scales = 'free_y') + 
  ylab( 'Parameter Estimate +/- 95% Bayesian Credible Intervals')

top_mod <- 
  top_mods %>% 
  select( vr, spp, climate_window, par , `50%`) %>% 
  spread( par, `50%`)

plot_2d_grid <- list()
obs_years <- list()

for( i in 1:nrow(top_mod)){ 
  vr <- top_mod$vr[i]
  spp <- top_mod$spp[i]
  window <- top_mod$climate_window[i]
  
  fn <- paste0('data/temp_data/', spp, '_growth_survival_dataframe.RDS')

  temp <- readRDS(fn)
  
  temp_mod <- top_mod[i,  ]
  
  if(is.na(window)){ 
    window <- 6
  }
  
  top_covs <- paste0(c('C.VWC.', 'C.T.'), window)
  
  all_years <- 
    temp %>%
    filter( Period == 'Historical') %>% 
    select(year, top_covs) %>% 
    distinct() 
  
  names( all_years ) <- c('year', 'moist', 'therm')
  
  vars <- expand.grid( temp_x = seq(-2, 3.5, by = 0.1), 
                       moist_y = seq(-2, 2, by = 0.1))

  plot_2d_grid[[i]] <- 
    data.frame(vars, temp_mod) %>% 
    mutate( response = temp_x*Temp + moist_y*Moist + temp_x*moist_y*Temp_x_Moist ) %>% 
    distinct( spp, vr, temp_x, moist_y, response )
  
  all_years$spp <- spp 
  all_years$vr <- vr 
  obs_years[[i]] <- all_years 
  
  if( all(is.na(plot_2d_grid[[i]]$response))){ 
    obs_years[[i]]$therm <- NA
    obs_years[[i]]$moist <- NA
  }
  
}


plot_2d_grid <- do.call(rbind, plot_2d_grid)
obs_years <- do.call(rbind, obs_years)

plot_2d_grid %>% group_by(spp, vr, moist_y, temp_x) %>% arrange( temp_x, moist_y)

gg2 <- 
  plot_2d_grid %>% 
  group_by(spp, vr) %>% 
  mutate( scaled_response = scale(response)) %>% 
  ggplot( aes( x = moist_y, y = temp_x, z = scaled_response, fill = scaled_response )) + 
  geom_tile() + 
  geom_contour(bins = 4) + 
  geom_point(data = obs_years, aes( x = moist, y = therm, z = NULL, fill = NULL ), alpha = 0.5) + 
  facet_grid(vr~spp) + 
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', name = 'response predicted') + 
  xlab( 'Annual Soil Moisture\nDry      --->       Wet') + 
  ylab( 'Annual Temperature\nCold      --->      Hot') + 
  ggtitle('Predicted effects of temperature and moisture in top models')

gg2 + coord_equal()

ggsave(gg2 + coord_equal(), 
       filename = 'figures/top_model_responses.png', width = 8, height = 6)

ggsave( 
  gg1 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + 
  xlab('parameter') + 
  ylab('Estimate +/- 95% Bayesian Credible Intervals') + 
  ggtitle( 'climate effects in top models'), 
  filename = 'figures/top_model_effects.png', width = 8, height = 6)

