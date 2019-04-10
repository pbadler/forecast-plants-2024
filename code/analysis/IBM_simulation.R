rm(list =ls())
library(tidyverse)

# functions ---------------------------------------------------------------------------- 

get_recruit_area <- function(spp) { 
  
  dat <- read.csv(paste0('lib/driversdata/data/idaho/speciesData/', as.character(spp), '/recSize.csv'))
  
  dat$area

}

#
load('code/figure_scripts/my_plotting_theme.Rdata')

quads <- read_csv('data/quad_info.csv')

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

fns <- dir('output/stan_fits', '.RDS', full.names = T)

fn_df <- 
  data.frame(fns) %>% 
  mutate( bname = basename(as.character(fns))) %>% 
  separate( bname, c('spp', 'vr', 'window', 'type'), sep = '_') %>% 
  mutate( type = str_remove(type, pattern = '\\.RDS$' )) 

fn_df <- 
  fn_df %>% 
  group_by( spp, vr , type ) %>% 
  arrange( window) %>% 
  mutate( rank = row_number() ) %>% 
  mutate( climate = rank == 1) %>% 
  select(-rank ) 

model_list <- 
  rbind( 
  fn_df %>% 
    filter( climate ) %>% 
    unite( type, vr, type, sep = '_') %>% 
    select(-window) %>% 
    spread( type, fns )
  , 
  fn_df %>% 
    select( - climate ) %>% 
    filter( window == 'none') %>% 
    unite( type , vr, type , sep = '_' ) %>% 
    select( -window) %>% 
    spread( type, fns) %>% 
    mutate( climate = F)
)


model_list <- 
  model_list %>% 
  group_by(spp) %>%  
  mutate( IBM_ID = row_number()) 

model_list %>% 
  write_csv('output/IBM_model_table.csv')

# loop species and climate / non-climate IBMs 
i = 1

for(i in 1:nrow( model_list)){ 
  spp <- model_list$spp[i]
  IBM_ID <- model_list$IBM_ID[i]
  
  gd <- readRDS( as.character(model_list$growth_data[i]))
  rd <- readRDS( as.character(model_list$recruitment_data[i]))
  sd <- readRDS( as.character(model_list$survival_data[i]))
  
  gfit <- readRDS( as.character(model_list$growth_model[i]))
  sfit <- readRDS( as.character(model_list$survival_model[i]))
  rfit <- readRDS( as.character(model_list$recruitment_model[i]))
  
  # generate predicted area per quad per year -------------- # 
  S <- binomial(link='logit')$linkinv(rstan::extract( sfit, 'IBM_mu')$IBM_mu)
  G <- rstan::extract( gfit, 'IBM_Y_hat')$IBM_Y_hat
  
  # Test that the predictions can be rescaled correctly to cm^2 #  
  Y_attrib <- gd$IBM_Y_attrib
  Y_center <- Y_attrib$`scaled:center`
  Y_scale <- Y_attrib$`scaled:scale`
  
  gdat <- readRDS(paste0('data/temp_data/', spp, '_growth_survival_dataframe.RDS'))
  
  Y1 <- (gd$IBM_Y*Y_scale + Y_center)
  Y2 <- gdat$logarea.t1
  all.equal(Y1, Y2)
  #------------------------------------------------------------ # 
  
  G <- exp( G*Y_scale + Y_center ) # tranform to cm scale
  R <- rstan::extract( rfit, 'IBM_Y_hat')$IBM_Y_hat
  
  a <- get_recruit_area(spp, my_path)
  
  K <- S*G  # survival by size 
  R <- R*median(a)
  
  K <- data.frame( quad = sd$IBM_quad_name, year = sd$IBM_year_name + 1, t(K))
  R <- data.frame( quad = rd$IBM_quad_name, year = rd$IBM_year_name + 1, t(R))
  
  K <- 
    K %>% 
    gather( sim, area, starts_with('X')) %>%
    group_by( quad, year, sim) %>% 
    summarise(area = sum( area ))
  
  R <- 
    R %>% 
    gather( sim, area, starts_with('X')) 
  
  A_pred <- 
    R %>% 
    left_join(K, by = c('quad', 'year', 'sim')) %>% 
    ungroup() %>% 
    gather( type, area, area.x, area.y) %>% 
    group_by( quad, year, sim) %>% 
    summarise( area = sum(area, na.rm = T)) %>%
    mutate( cover = 100*area/(100*100))
  
  A_pred <- 
    expand.grid( year = seq( min( A_pred$year) - 1, max(A_pred$year) ) + 1, quad = unique( A_pred$quad )) %>% 
    left_join(A_pred, by = c('quad', 'year')) %>% 
    left_join(quads, by = c('quad' = 'QuadName')) %>% 
    select( quad, Treatment, year, cover, sim )
  
  A_pred$spp <- spp 
  
  saveRDS(A_pred, paste0( 'output/IBM_', spp, '_model_', IBM_ID, '_cover_predictions.RDS'))
}  
