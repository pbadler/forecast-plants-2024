rm(list = ls())
library(tidyverse)
library(climwin)
library(scales)
library(lubridate)

c_trans <- function(a, b, breaks = b$breaks, format = b$format) {
  a <- as.trans(a)
  b <- as.trans(b)
  name <- paste(a$name, b$name, sep = "-")
  trans <- function(x) a$trans(b$trans(x))
  inv <- function(x) b$inverse(a$inverse(x))
  
  trans_new(name, trans, inverse = inv, breaks = breaks, format=format)
  
}

label_func <- function(x) { 
  paste( month.abb[month(x)], mday(x), sep = '-')
}

load_to_name <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

rev_date <- c_trans("reverse", "date")

results_files <- dir('data/temp_data/', pattern = 'weekly_ClimWin.*rda', full.names = T)

results_files <- data.frame( file_name = results_files ) %>% 
    mutate( species = str_extract( file_name , '[A-Z]{4}'), 
            vital_rate = str_extract(file_name, '(growth)|(survival)'), 
            VWC_cov = str_detect( file_name, 'VWC_cov')) %>%
  filter( !is.na(vital_rate), 
          !(vital_rate == 'survival' & !VWC_cov & !str_detect(file_name, 'glmer'))) %>% 
  arrange(species, vital_rate, VWC_cov)


window_labels <- 
  data.frame( date = as_date( ymd('1997-01-01'):ymd('1999-06-15')) ) %>% 
  mutate( mday = mday(date),
          week = week(date), 
          month = month(date), 
          year_abs = year(date) - min(year(date))) %>% 
  mutate( week = ifelse(week == 53, 52, week)) %>% 
  mutate( week_abs = week + year_abs*52 ) %>% 
  mutate( weeks_back = max(week_abs) - week_abs, 
          years_back = max(year_abs) - year_abs) %>% 
  mutate( window_label = paste( month.abb[month], mday, years_back, 'year back')) %>%
  mutate( window_label = str_remove(window_label, '0 year back')) %>% 
  group_by( year_abs, week_abs) %>% 
  filter( date == min(date)) %>% 
  ungroup() %>% 
  select( date, weeks_back, window_label)


for( i in 1:nrow( results_files )){
  # Choose species and vital rate 
  
  temp_species <- results_files$species[i]
  temp_vital_rate <- results_files$vital_rate[i]
  temp_with_VWC <- results_files$VWC_cov[i]
  res <- load_to_name(results_files$file_name[i])

  for( j in 1:nrow(res$combos)){ 
    temp_stat <- res$combos$stat[j]
    temp_func <- res$combos$func[j]
    temp_var <- res$combos$climate[j]
    temp_dat <- res[[j]]$Dataset
    
    
    temp_dat <- 
      temp_dat %>% 
      left_join(window_labels, by = c('WindowOpen' = 'weeks_back')) %>%
      rename( WindowOpenLabel = window_label, 
              foo_dateOpen = date) %>%
      left_join(window_labels, by = c('WindowClose' = 'weeks_back')) %>% 
      rename( WindowCloseLabel = window_label, 
              foo_dateClose = date)
    
    
    climWin_plot <- plotdelta(temp_dat)
    plot_info <- ggplot_build( plotdelta(temp_dat) )
    climWin_theme <- plot_info$plot$theme
    
    best_window <- temp_dat[1, ] %>% 
      mutate( label = paste( WindowOpenLabel, WindowCloseLabel, sep = ' to ')) %>% 
      mutate( label = str_remove_all(label, "[0-9]+ year back")) %>%
      mutate( label = paste(label, round(ModelBeta, 2), sep = ': Beta = '))
    
    midAIC <- sum(range(temp_dat$deltaAICc))/2
    midBeta <- sum(range(temp_dat$ModelBeta))/2
    
    temp_title <- paste( temp_species, temp_vital_rate, temp_func, temp_stat, temp_var, sep = ' ')
    if( temp_with_VWC ) { 
      temp_title <- paste( temp_title, '\n (with best VWC)')
    }
    
    aic_plot <- 
      temp_dat %>% 
        ggplot( aes( y = foo_dateOpen, x = foo_dateClose, fill = deltaAICc)) + 
        geom_tile() + 
        geom_point(data = best_window) + 
        ggrepel::geom_label_repel(data = best_window, aes(label = label, fill = NULL), min.segment.length = Inf,
                                 point.padding = 1, box.padding = 1, nudge_x = 50, nudge_y = -20) + 
        scale_fill_gradient2(high = 'blue', low = 'red', mid = 'yellow', midpoint = midAIC ) +
        scale_x_continuous(trans = rev_date, labels = label_func, name = 'Window Close') + 
        scale_y_continuous(trans = rev_date, labels = label_func, name = 'Window Open') + 
        climWin_theme  + 
        ggtitle(temp_title)

    beta_plot <- 
      temp_dat %>% 
        ggplot( aes( y = foo_dateOpen, x = foo_dateClose, fill = ModelBeta)) + 
        geom_tile() + 
        geom_point(data = best_window) + 
        ggrepel::geom_label_repel(data = best_window, aes(label = label, fill = NULL), min.segment.length = Inf,
                                  point.padding = 1, box.padding = 1, nudge_x = 50, nudge_y = -20) + 
        scale_fill_gradient2(high = 'blue', low = 'red', mid = 'yellow', midpoint = midBeta ) +
        scale_x_continuous(trans = rev_date, labels = label_func, name = 'Window Close') + 
        scale_y_continuous(trans = rev_date, labels = label_func, name = 'Window Open') + 
        climWin_theme  + 
        ggtitle(temp_title)
    
    temp_title <- str_replace_all( str_squish( temp_title), c(' ' = '_', '\\('='', '\\)'=''))
    aic_output_figure <- paste0( 'data/temp_data/', 'AIC_figure_', temp_title, '.rds')
    beta_output_figure <- paste0( 'data/temp_data/', 'beta_figure_', temp_title, '.rds')
    
    saveRDS(aic_plot, file = aic_output_figure)
    saveRDS(beta_plot, file = beta_output_figure)
  }
}

aic_figs <- dir( '~/Dropbox/projects/forecast-plants/data/temp_data/temp_figure_files', pattern = 'AIC_figure_', full.names = T)
beta_figs <- dir( '~/Dropbox/projects/forecast-plants/data/temp_data/temp_figure_files', pattern = 'beta_figure_', full.names = T)

all_results <- 
  data.frame( fig_file = c( aic_figs, beta_figs) ) %>% 
  mutate( species = str_extract(fig_file , '[A-Z]{4}')) %>% 
  mutate( type = str_extract( fig_file , '(AIC|beta)')) %>% 
  mutate( vr = str_extract( fig_file, '(growth|survives)')) %>% 
  mutate( log = str_detect(fig_file, 'log')) %>%
  mutate( stat = str_extract( fig_file, 'sum|mean')) %>% 
  mutate( clim = str_extract(fig_file, '(VWC)|(GDD)|(PRCP)|(TMAX)|(TAVG)|(TAVG_0)|(TDIF)')) %>% 
  mutate( with_VWC_covariate = str_detect( fig_file, 'with_best_VWC'))

saveRDS(all_results, 'data/temp_data/all_results.rds')
file.copy( 'data/temp_data/all_results.rds', to = 'compare_windows/data', overwrite = T)
file.copy(from = all_results$fig_file, to = 'compare_windows/data', overwrite = T)
