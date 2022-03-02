rm(list = ls())

library(tidyverse)

load('output/growth_model_combos.rda')
my_path <- 'output/model_scores/'
vr <- 'growth'

get_scores <- function(my_path, vr ){ 
  
  files <- dir(path = my_path, 
               pattern = paste0( '^', vr, '*'), 
               full.names = T)
  
  scores <- lapply( files, read_rds )
  
  scores <- 
    lapply( scores, function(x){ 
    x$min_rhat = min(x$rhat); 
    x$max_rhat = max(x$rhat); 
    return(x[ -( which(names(x) == 'rhat'))]) } )

  scores <- do.call(bind_rows, scores)
  scores$vr <- vr 
  return( scores )
}  

scores <- get_scores('output/model_scores/', 'growth')

scores %>% 
  group_by( species, model ) %>% 
  summarise( ins_lppd = sum( ins_lppd ), 
             oos_lppd = sum( oos_lpd ), 
             oos_sse = sum(oos_sse)) %>% 
  mutate( oos_lppd =  oos_lppd/max(oos_lppd)) %>% 
  mutate( oos_sse =  oos_sse/min(oos_sse )) %>% 
  gather( stat, value, ins_lppd:oos_sse ) %>% 
  filter( stat != "ins_lppd") %>% 
  ggplot( aes( x = model, y = value, group = paste( species, stat), color = stat)) + 
  geom_line() + 
  facet_wrap( ~ species, scale = 'free_x' ) + 
  coord_flip() + 
  theme_bw() + 
  ylab('out of sample support')


