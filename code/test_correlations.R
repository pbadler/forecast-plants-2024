rm(list = ls())

library(tidyverse)
library(imputeTS)
library(astsa)
library(forecast)
library(lme4)
library(zoo)
library(lubridate)
library(climwin)

my_rle_fun <- function(x) { unlist( sapply( rle( x)$lengths, function(x) rep(x, x)))  } 

get_gaps <- function(x){ 
  gap_vals <- NA
  for( i in 1:length(x) ){ 
    gap_vals[i] <- all( is.na(x[i]) , any( !is.na( x[1:(i-1)] )) , any( !is.na( x[(i+1):length(x)] ) ))
  }
  return(gap_vals)
}

load( file = 'data/temp_data/st_anom.rda')

sp_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

species <- 'POSE'

intra <- paste0( 'W.' , species )
dat <- readRDS( paste0 ( 'data/temp_data/', paste0( species, '_size.RDS')))
dat <- dat[ dat$year < 1960 , ]

mnyear <- min(dat$year)
mxyear <- max(dat$year)
pid <- unique( dat$pid )

dat <- 
  expand.grid( year = mnyear:mxyear, pid = unique(pid)) %>% 
  left_join(dat, by = c('pid', 'year'))  %>% 
  arrange(pid, year) 

my_dat <- 
  dat %>% 
  group_by( pid ) %>% 
  mutate( age = ifelse( year == min(year[!is.na(area)]), 1, age )) %>% 
  dplyr::select( pid, age, year, area , starts_with('W.'))

# find short sequences of missing values
# interpolate these
# fit time series models to each plant with long sequence
# look at residuals 
# look for year effects in the residuals 

good_dat <- 
  my_dat %>% 
  group_by( pid ) %>% 
  arrange( pid, year ) %>% 
  mutate( gap = get_gaps( area ) ) %>% 
  mutate( gap_length = ifelse( gap, my_rle_fun(gap) , 0 )) %>% 
  filter( !is.na(area) | ( gap & gap_length < 10)) %>% 
  mutate( gap = as.numeric( gap )) %>% 
  select( - gap_length) %>% 
  mutate( age = row_number()) %>%
  group_by( pid ) %>%   
  arrange( pid, year ) %>% 
  mutate( missing = is.na(area)) %>% 
  mutate( area = log(area)) %>% 
  mutate_at( vars('area', starts_with('W')), ts) %>% 
  mutate_at( vars('area', starts_with('W')), na_interpolation) %>%
  select( year, pid, age, area, starts_with('W'), missing) %>% 
  group_by( pid ) %>% 
  mutate( area0 = lag(area))


long_series <- 
  good_dat %>% 
  group_by( pid ) %>% 
  mutate( n_years = max(year) - min(year)) %>% 
  filter( n_years > 25 ) %>% 
  mutate( plot = str_extract(pid, 'Q\\d+')) %>% 
  group_by( pid ) %>% 
  mutate( growth = area - area0) %>% 
  mutate( size2 = scale( area )) #%>% 
  #filter( plot %in% c('Q19', 'Q20', 'Q7', 'Q8'))

long_series %>% 
  ggplot( aes( x = year, y = area, group = pid )) + 
  geom_line() + 
  facet_wrap( ~ plot )

long_series %>% 
  ggplot( aes( x = year, y = size2, group = pid )) + 
  geom_line() + 
  facet_wrap( ~ plot ) + 
  geom_point(data = long_series %>% filter( missing ) , color = 'red' )

long_series %>% 
  group_by( pid) %>% 
  summarise( va = var( area), ma = mean( area),  mg = mean(growth, na.rm = T),  vg = var( growth, na.rm = T) , mc = mean( W.PSSP + W.POSE + W.HECO)) %>% 
  ggplot( aes( x = ma , y = mg)) + 
  geom_point() 

long_series %>% 
  ggplot( aes( x = year, y = growth, group = pid )) + 
  geom_line() + 
  stat_summary( fun = 'mean', aes( group = 1), geom = 'line', color = 'red') 

library(lme4)

m1 <- lmer( data = long_series, area ~ area0 + eval(parse( text = intra)) + (1|pid) + (1|year))

mg <- lmer( data = long_series, growth ~ eval(parse(text = intra)) + (1|year) + (1|plot)  )
summary( mg )

plot( 1924:1957, ranef(m1)$year[,1] , type = 'l')
points( 1924:1957, ranef(mg)$year[,1], type = 'l', col = 'red', lty = 2)

last_summer_clim <- 
  st_anom %>% 
  filter( end == 8 , len == 6) %>% 
  ungroup() %>% 
  select( var, year, rllm )

winter_clim <- 
  st_anom %>% 
  filter( end == 1, len == 6) %>% 
  ungroup() %>% 
  select( var, year , rllm )

spring_clim <- 
  st_anom %>% 
  filter( end == 0, len == 3) %>% 
  ungroup() %>% 
  select( var, year, rllm)

growth <- 
  last_summer_clim %>% 
  group_by( var) %>% 
  mutate( E_S = scale(rllm)) %>%
  left_join(winter_clim, by = c('year', 'var')) %>% 
  mutate( E_W = scale(rllm.y)) %>% 
  left_join( spring_clim, by = c('year', 'var')) %>% 
  mutate( E_Sp = scale( rllm )) %>% 
  left_join( 
    (data.frame( ranef(m1)$year ) %>% mutate( year = as.numeric( row.names(.))) %>% rename( g1 = X.Intercept.) %>% select( year, g1 ) )) %>% 
  left_join(
    (data.frame( ranef(mg)$year ) %>% mutate( year = as.numeric( row.names(.))) %>% rename( g2 = X.Intercept.) %>% select( year, g2 ) )
  ) %>% 
  mutate( g1 = scale( g1 ), g2 = scale(g2)) %>% 
  gather( model, value , g1:g2)

growth %>% 
  ggplot( aes( x = E_W, y = value, color = var)) + 
  geom_point() + 
  geom_smooth(se = F, method = 'lm') +
  scale_color_manual(values = c('blue', 'red', 'purple')) + 
  facet_wrap(~ model )

growth %>% 
  ggplot( aes( x = E_S, y = value, color = var)) + 
  geom_point() + 
  geom_smooth(se = F, method = 'lm') +
  scale_color_manual(values = c('blue', 'red', 'purple')) + 
  facet_wrap(~ model )

growth %>% 
  ggplot( aes( x = E_Sp, y = value, color = var)) + 
  geom_point() + 
  geom_smooth(se = F, method = 'lm') +
  scale_color_manual(values = c('blue', 'red', 'purple')) + 
  facet_wrap(~ model )

growth %>% 
  filter( var == 'PRCP') %>% 
  ggplot( aes( x = year,  y = value, linetype = model )) + 
  geom_line() + 
  geom_line(aes( y = E_W, color = var ), linetype = 1, size = 1) + 
  scale_color_manual(values = c('blue', 'red', 'purple'))

growth %>% 
  filter( var == 'PRCP') %>% 
  ggplot( aes( x = year,  y = value, linetype = model )) + 
  geom_line() + 
  geom_line(aes( y = E_S, color = var ), linetype = 1, size = 1) + 
  scale_color_manual(values = c('blue', 'red', 'purple'))

growth %>% 
  filter( var == 'TAVG') %>% 
  ggplot( aes( x = year,  y = value, linetype = model )) + 
  geom_line() + 
  geom_line(aes( y = E_S, color = var ), linetype = 1, size = 1) + 
  scale_color_manual(values = c('blue', 'red', 'purple'))

# using my climate vars : 

mod_with_intra <- lmer( data = good_dat, area ~ area0 + eval( parse(text = intra) ) + ( 1 | year))
summary(mod_with_intra)
ranef(mod_with_intra)

plot( min(good_dat$year):(max(good_dat$year)-1), ranef(mod_with_intra)$year$`(Intercept)`, type = 'l' ,
      xlab = 'year', ylab = 'random intercept')
points( min(good_dat$year):(max(good_dat$year)-1), ranef(mod_with_intra)$year$`area0`, type = 'l' ,
      xlab = 'year', ylab = 'random intercept', lty = 2)
legend('topright', c('Year', 'Year x area') , lty = c(1,2), title = 'Year Effect')

mod_without_intra <- lmer( data = good_dat, area ~ area0 + (area0|year))
summary(mod_without_intra)
ranef(mod_without_intra)

plot( min(good_dat$year):(max(good_dat$year)-1), ranef(mod_with_intra)$year$`(Intercept)`, type = 'l' ,
      xlab = 'year', ylab = 'random intercept')
points( min(good_dat$year):(max(good_dat$year)-1), ranef(mod_with_intra)$year$`area0`, type = 'l' ,
        xlab = 'year', ylab = 'random intercept', lty = 2)
points( min(good_dat$year):(max(good_dat$year)-1), ranef(mod_without_intra)$year$`(Intercept)`, type = 'l' ,
      xlab = 'year', ylab = 'random intercept', col = 'red')
points( min(good_dat$year):(max(good_dat$year)-1), ranef(mod_without_intra)$year$`area0`, type = 'l' ,
        xlab = 'year', ylab = 'random intercept', lty = 2, col = 'red')

legend('topright', c('Year', 'Year x area', 'Year w/o comp.', 'Year x area w/o comp.') , lty = c(1,2, 1,2), col = c(1,1, 2,2), title = 'Year Effect')

ranef_df <- data.frame( ranef(mod_with_intra)$year)
ranef_df$year <- as.numeric (row.names(ranef_df) )
ranef_df$ranef <- ranef_df$X.Intercept.

year_weights <- 
  good_dat %>% 
  group_by( year ) %>%
  summarise( n = n())

test_correlation_df <- 
  ranef_df %>% 
  left_join(st_anom) %>% 
  filter( !is.na(var) ) %>%  
  ungroup() %>% 
  filter( end <= 12, len <= 12) %>% 
  group_by(var, end, len) %>%
  arrange( var, year ) %>% 
  filter( !is.na(rllm), !is.na(ranef)) 

cor_results <- 
  test_correlation_df %>% 
  summarise( r = cor.test( ranef , rllm )$estimate, n_years = n() )

label_best <-
  cor_results %>% 
  mutate( dummy_date = ymd( paste( 1940, 6, 1)) , 
          end_month = month.name[ month( dummy_date - months(end) ) ],
          start_month = month.name[ month( dummy_date - months(end + len))], 
          years_back = end %/% 12 ) %>% 
  group_by( var ) %>% 
  filter( r > quantile(r, 0.95) | r < quantile(r, 0.05)) %>%
  mutate( start_years_back = (end + len) %/% 12, end_years_back = (end) %/% 12 ) %>% 
  mutate( label = paste( paste( start_month, start_years_back, 'years back to', end_month, end_years_back, 'years back'))) %>% 
  mutate( r = round( r, 2)) %>% 
  ungroup() %>% 
  mutate( letter = LETTERS[row_number()] )



cor_results %>% 
  ggplot( aes( x = end, y = len, z = r, fill = r)) + 
  geom_raster(interpolate = F) + 
  facet_wrap(~ var) + 
  scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue', limits = c(-1,1)) + 
  geom_point(data = label_best, aes( x = end, y = len), shape = 19, color = 'black') + 
  ggrepel::geom_label_repel(data = label_best, aes( x = end, y = len, label = paste0('r = ', r)),  force = 2, hjust = 1) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  xlab( 'Window End (Months Back)' )  + 
  ylab( 'Window Duration (Months)')

label_best %>% View

ggsave(filename = paste0('figures/', species, 'clim_correlation.pdf'), 
       cor_results %>% 
        ggplot( aes( x = end, y = len, z = r, fill = r)) + 
        geom_raster(interpolate = T) + 
        geom_contour(bins = 3, show.legend = T) + 
        facet_wrap(~ var) + 
        scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue') + 
        geom_point(data = label_best, aes( x = end, y = len), shape = 1, color = 'white') + 
        ggrepel::geom_text_repel(data = label_best, aes( x = end, y = len, label = paste0(letter, ': r = ', r)), color = 'white',  force = 2, hjust = 1) + 
        theme_bw() + 
        theme(panel.grid = element_blank()) + 
        xlab( 'Window End (Months Back)' )  + 
        ylab( 'Window Duration (Months)')
)


label_best %>% 
  select( var, letter, r, label ) 

cor_results %>% 
  group_by( var ) %>% 
  filter( abs(r) == max(abs(r))) 

test_correlation_df %>% 
  left_join(cor_results, by = c('var', 'end', 'len', 'var')) %>% 
  group_by(var) %>% 
  filter( abs(r) == max(abs(r))) %>% 
  left_join(year_weights) %>% 
  ggplot( aes( x = rllm,  y = ranef , size = n)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth( method = 'lm', se = F) + 
  ggrepel::geom_text_repel( aes(label = year), size = 3, alpha = 0.5)  + 
  facet_wrap(~ var, scales = 'free' ) + 
  ylab( 'Growth Year Effect' ) + 
  xlab( 'Climate Value') + 
  scale_size_continuous(range = c(0.5,4), guide = 'none')

top_covariates <- 
  cor_results %>% 
  group_by( var ) %>% 
  filter( abs(r) == max(abs(r))) 

top_covariates

fit_dat <- 
  st_anom %>% 
  left_join(top_covariates) %>% 
  group_by( var ) %>% 
  filter( !is.na(r)) %>%
  ungroup() %>% 
  select( year, var, rllm) %>% 
  spread( var, rllm )  %>% 
  left_join(good_dat, by = 'year') %>% 
  filter( !is.na(PRCP), !is.na(TAVG), !is.na(VWC))


m1 <- lmer(data = fit_dat, 
           area ~ area0 + eval(parse(text = intra)) + PRCP + (area0|year))

summary( m1 )

ranef_m1 <- data.frame( ranef(m1)$year )
ranef_m1$year <- row.names(ranef_m1)
ranef_m0 <- data.frame( ranef(mod_with_intra)$year )
ranef_m0$year <- row.names( ranef_m0)

var( ranef_m1$X.Intercept. )
var( ranef_m0$X.Intercept.)

plot( ranef_m0$year, ranef_m0$X.Intercept. , type = 'l', xlab = 'Year', ylab = 'Year Effect on Growth')
points( ranef_m1$year, ranef_m1$X.Intercept., type = 'l', col = 'red')
legend('topright', c('no climate', 'with climate') , lty = 1, col = c(1,2), title = 'Model' )

