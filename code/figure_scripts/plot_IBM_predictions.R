rm(list = ls())

IBM_table <- read_csv('output/IBM_model_table.csv')

cover_obs <- readRDS('data/temp_data/all_cover.RDS') # observed cover per quadrat 

cover_obs <- 
  cover_obs %>% 
  rename('quad' = QuadName, 'cover_obs' = cover) %>% 
  ungroup() %>% 
  mutate( era = cut(year, include.lowest = T, breaks = c(min(year)-1, 1960, 2004, 2018), 
                    labels = c('early', 'mid', 'late'))) 
  
era_labels <- 
  cover_obs %>% 
  distinct(year, era)

sppList <- unique(cover_obs$spp)

i <- 1 

spp <- IBM_table$spp[i] 
IBM_ID <- IBM_table$IBM_ID[i] 

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


q_pred_df 

t_cover_obs <- 
  cover_obs %>% 
  ungroup %>%
  filter( spp == !! spp, era != 'mid') %>% 
  group_by( year, Treatment, spp, era )  %>% 
  summarise( cover_obs = mean(cover_obs, na.rm =T)) 

t_cover_pred %>% 
  filter( era != 'mid') %>% 
  ggplot( aes( x = year, y = med50, color = Treatment)) + 
  geom_ribbon(aes( ymin = low25, ymax = upper75, fill = Treatment), alpha = 0.2) + 
  geom_line() + 
  geom_point(data = t_cover_obs, aes(y = cover_obs)) + 
  facet_grid( Treatment ~ era, scale = 'free_x')


