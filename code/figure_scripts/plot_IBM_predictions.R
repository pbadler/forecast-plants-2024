rm(list = ls())

library(tidyverse)
load('code/figure_scripts/my_plotting_theme.Rdata')

IBM_table <- read_csv('output/IBM_model_table.csv')
cover_obs <- readRDS('data/temp_data/all_cover.RDS') # observed cover per quadrat 

# x-axis labeling function --------------------------- # 
my_breaks = function(x) { 
  if (x < 2006) { 
    seq(ceiling(1925), floor(x[2]), by = 10) 
  }else if ( x > 2006) { 
    seq(ceiling(x[1]), floor(2016), by = 2) 
  }
}
# ---------------------------- # 

cover_obs <- 
  cover_obs %>% 
  rename('quad' = QuadName, 'cover_obs' = cover) %>% 
  ungroup() %>% 
  mutate( era = cut(year, include.lowest = T, breaks = c(min(year)-1, 1961, 2006, 2018), 
                    labels = c('early', 'mid', 'late'))) 
  
era_labels <- 
  cover_obs %>% 
  distinct(year, era)

sppList <- unique(cover_obs$spp)


for( i in 1:nrow(IBM_table)){ 

  spp <- IBM_table$spp[i] 
  IBM_ID <- IBM_table$IBM_ID[i] 
  
  type <- ifelse( IBM_table$climate[i], 'climate', 'none')
  
  cover_pred <- readRDS( paste0( 'output/IBM_', spp, '_model_', IBM_ID, '_cover_predictions.RDS') )
  
  q_cover_pred <- 
    cover_pred %>% 
    group_by( year, quad, Treatment, spp) %>% 
    summarise( avg = mean(cover, na.rm = T), 
               low5 = quantile(cover, 0.05, na.rm = T), 
               low25 = quantile(cover, 0.25, na.rm = T),
               med50 = quantile(cover, 0.5, na.rm = T),
               upper75 = quantile(cover, 0.75, na.rm = T), 
               upper95 = quantile(cover, 0.95, na.rm = T)) %>% 
    mutate( year_label = as.numeric( str_sub(year, 3, 5)) ) %>%
    ungroup() %>% 
    left_join(era_labels, by = 'year')
  
  t_cover_pred <- 
    cover_pred %>% 
    group_by( year, sim, Treatment, spp ) %>% 
    summarise( cover = mean(cover, na.rm = T)) %>% 
    group_by( year, Treatment, spp) %>% 
    summarise( avg = mean(cover, na.rm = T), 
               low5 = quantile(cover, 0.05, na.rm = T), 
               low25 = quantile(cover, 0.25, na.rm = T),
               med50 = quantile(cover, 0.5, na.rm = T),
               upper75 = quantile(cover, 0.75, na.rm = T), 
               upper95 = quantile(cover, 0.95, na.rm = T)) %>% 
    mutate( year_label = as.numeric( str_sub(year, 3, 5)) ) %>% 
    ungroup() %>% 
    left_join(era_labels, by = 'year')
  
  q_pred_df <- 
    q_cover_pred %>% 
    left_join( cover_obs, by = c('quad', 'year', 'Treatment', 'spp', 'era') )
  
  t_cover_obs <- 
    cover_obs %>% 
    ungroup %>%
    filter( spp == !! spp, era != 'mid') %>% 
    group_by( year, Treatment, spp, era )  %>% 
    summarise( cover_obs = mean(cover_obs, na.rm =T)) 
  
  
  gg_cover_pred <- 
    t_cover_pred %>% 
    filter( era != 'mid') %>% 
    ggplot( aes( x = year, y = med50, color = Treatment)) + 
    geom_ribbon(aes( ymin = low5, ymax = upper95, fill = Treatment), alpha = 0.2, color = NA) + 
    geom_ribbon(aes( ymin = low25, ymax = upper75, fill = Treatment), alpha = 0.4, color = NA) + 
    geom_line() + 
    geom_point(data = t_cover_obs, aes(y = cover_obs)) + 
    facet_grid( Treatment ~ era, scale = 'free') + 
    ylab('Cover predicted/observed (%)') + 
    xlab('Year') + 
    ggtitle(paste0( IBM_table$spp[i],' ', type, ' model' )) + 
    scale_fill_manual(values = my_colors[2:4]) + 
    scale_color_manual(values = my_colors[2:4]) + 
    scale_x_continuous(breaks = my_breaks) + 
    my_theme 

  saveRDS(gg_cover_pred, paste0( 'output/gg_cover_pred_', IBM_table$spp[i], '_', type, '.RDS'))
  
}  
