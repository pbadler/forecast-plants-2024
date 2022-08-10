rm(list = ls())

library(climwin)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(ggpubr)
library(lme4)

source('code/analysis/functions.R')

result_files <- dir('output/growth_models', pattern = 'growth.*mer_monthly_ClimWin.rda', full.names = T)

for( i in result_files){ 
  load(i)
}

growth_windows <- bind_rows( 
  makeWindowTable(ARTR_growth_monthly_ClimWin, "ARTR", "growth"), 
  makeWindowTable(ARTR_growth_no_intxn_monthly_ClimWin, 'ARTR', 'growth_no_intxn'), 
  makeWindowTable(ARTR_growth_small_monthly_ClimWin, 'ARTR', 'growth_small'), 
  
  makeWindowTable(HECO_growth_monthly_ClimWin, "HECO", "growth"), 
  makeWindowTable(HECO_growth_no_intxn_monthly_ClimWin, 'HECO', 'growth_no_intxn'), 
  makeWindowTable(HECO_growth_small_monthly_ClimWin, 'HECO', 'growth_small'), 

  makeWindowTable(POSE_growth_monthly_ClimWin, "POSE", "growth"), 
  makeWindowTable(POSE_growth_no_intxn_monthly_ClimWin, 'POSE', 'growth_no_intxn'), 
  makeWindowTable(POSE_growth_small_monthly_ClimWin, 'POSE', 'growth_small'), 
  
  makeWindowTable(PSSP_growth_monthly_ClimWin, "PSSP", "growth"), 
  makeWindowTable(PSSP_growth_no_intxn_monthly_ClimWin, 'PSSP', 'growth_no_intxn'), 
  makeWindowTable(PSSP_growth_small_monthly_ClimWin, 'PSSP', 'growth_small'), 
  
  
) 
write_csv( growth_windows, 'output/growth_models/top_growth_windows.csv')

result_files_s <- dir('output/survival_models', pattern = 'survival.*mer_monthly_ClimWin.rda', full.names = T)

for( i in result_files_s){ 
  load(i)
}

survival_windows <- bind_rows(
  makeWindowTable(ARTR_survival_monthly_ClimWin, 'ARTR', 'survival'), 
  makeWindowTable(ARTR_survival_no_intxn_monthly_ClimWin, 'ARTR', 'survival_no_intxn'), 
  
  makeWindowTable(HECO_survival_monthly_ClimWin, 'HECO', 'survival'), 
  makeWindowTable(HECO_survival_no_intxn_monthly_ClimWin, 'HECO', 'survival_no_intxn'), 
  
  makeWindowTable(POSE_survival_monthly_ClimWin, 'POSE', 'survival'), 
  makeWindowTable(POSE_survival_no_intxn_monthly_ClimWin, 'POSE', 'survival_no_intxn'), 
  
  makeWindowTable(PSSP_survival_monthly_ClimWin, 'PSSP', 'survival'), 
  makeWindowTable(PSSP_survival_no_intxn_monthly_ClimWin, 'PSSP', 'survival_no_intxn')
)

write_csv( survival_windows, 'output/survival_models/top_survival_windows.csv')

spList <- c('ARTR', 'HECO', 'POSE', 'PSSP')

growth_windows <- 
  growth_windows %>% 
  split( .$type )

survival_windows <- 
  survival_windows %>% 
    split( .$type) 

# Make plots for each species
for( sp in spList){ 
  make_window_plots(growth_windows$growth, sp )
  make_window_plots(survival_windows$survival, sp )
}


# Plots without size x climate interaction  
for( sp in spList){ 
  plot_windows_no_intxn(growth_windows$growth_no_intxn, sp)
  plot_windows_no_intxn(growth_windows$growth_small, sp )
  plot_windows_no_intxn(survival_windows$survival_no_intxn, sp )
}
