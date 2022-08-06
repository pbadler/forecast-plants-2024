rm(list = ls() )

# 
# Produce the dataframes with climate and demographic data for stan 
# If the temporary demographic data files are missing then the "get_all_demographic_data" script needs to be re-run 
#

library(tidyverse)

source( 'code/prepare_data/prep_vital_rate_df_functions.R')


prep_vital_rate_df <- function( clim_file = 'data/temp/all_clim_covs.csv', gfile , sfile , rfile = NULL){ 
  
  clim_vars <- scan(file = file.path('data', 'select_clim_vars.txt'), what = 'char')
  clim <- scale_climate_vars(clim_file, clim_vars)
  
  if(!(gfile == '') & !(sfile == '')){ 
    df <- clean_growth_survival(clim, gfile, sfile)
  }else if(!(rfile == '')){
    df <- clean_recruitment(clim, rfile)
  }else{
    stop('too many files provided')
  }
  
  fname <- paste( unique(unlist(
    lapply(c(gfile, sfile, rfile), 
           function(x)
           {
             strsplit(sub('\\..*$', '', basename(x), perl= T), '_')
           }
    )
  )), collapse = '_')
  
  write_csv(df, file = file.path( 'data', 'temp', paste0(fname, '_dataframe.csv')))
  
} 

data_df <- expand.grid( species = c('ARTR', 'HECO', 'POSE', 'PSSP'), 
                        vr = c('growth_survival', 'recruitment'))


data_df <- 
  data_df %>% 
  mutate( gfile = ifelse( vr == 'growth_survival', 
                          paste0('data/temp/', species, '_growth.csv'), 
                          ''), 
          sfile = ifelse( vr == 'growth_survival',
                          paste0('data/temp/', species, '_survival.csv'),
                          ''), 
          rfile = ifelse( vr == 'recruitment',
                          paste0('data/temp/', species, '_recruitment.csv'),
                          ''))
i <- 1
for( i in 1:nrow(data_df)){
  data_df$gfile[i]
  
  prep_vital_rate_df(gfile = data_df$gfile[i], sfile = data_df$sfile[i], rfile = data_df$rfile[i])
  
}

