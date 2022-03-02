rm(list = ls())

library(tidyverse)

score_fls <- dir('output', pattern = 'model_scores.RDS', full.names = T)

all_ranks <- lapply(score_fls, readRDS ) 

rank_mods <- do.call( rbind, all_ranks)

rank_mods %>% 
  write_csv('output/model_ranks.csv')


# plots -------------------------------------------------------------------- # 

rank_mods %>% 
  mutate(window_lab = ifelse (is.na(climate_window), 0, climate_window)) %>% 
  group_by( spp, vr) %>% 
  mutate( top = oos_lppd == max(oos_lppd)) %>% 
  ungroup() %>% 
  ggplot( aes( x = climate_window, y = oos_lppd, color = spp)) + 
  geom_line() + 
  geom_point(data = . %>% filter(top), aes( x = climate_window, y = oos_lppd), shape = 2, show.legend = F) + 
  facet_wrap( vr ~ spp, scales = 'free_x') + 
  coord_flip()

rank_mods %>% 
  mutate(window_lab = ifelse (is.na(climate_window), 0, climate_window)) %>% 
  group_by( spp, vr) %>% 
  mutate( top = oos_mse == min(oos_mse)) %>% 
  ggplot( aes( x = window_lab , y = oos_mse, color = spp)) + 
  geom_line() + 
  geom_point(data = . %>% filter(top), aes( x = window_lab, y = oos_mse), shape = 2, show.legend = F) + 
  facet_wrap( vr ~ spp, scales = 'free_x') + 
  coord_flip()

rank_mods %>% 
  ggplot( aes ( x = oos_lppd, y = oos_mse)) + 
  geom_point() + 
  geom_smooth(se = F, method = 'lm' ) + 
  facet_wrap(vr ~ spp, scales = 'free' )

