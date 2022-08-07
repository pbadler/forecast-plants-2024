rm(list = ls())

library(climwin)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(ggpubr)
library(lme4)

source('code/analysis/functions.R')

result_files <- dir('output/growth_models/', pattern = 'growth_mer_monthly_ClimWin.rda', full.names = T)

for( i in result_files){ 
  load(i)
}

growth_windows <- bind_rows( 
  makeWindowTable(ARTR_growth_monthly_ClimWin, "ARTR", "growth"), 
  makeWindowTable(HECO_growth_monthly_ClimWin, 'HECO', 'growth'), 
  makeWindowTable(POSE_growth_monthly_ClimWin, 'POSE', 'growth'), 
  makeWindowTable(PSSP_growth_monthly_ClimWin, 'PSSP', 'growth')
) 

write_csv( growth_windows, 'output/growth_models/top_growth_windows.csv')

result_files_s <- dir('output/survival_models/', pattern = 'survival_mer_monthly_ClimWin.rda', full.names = T)
for( i in result_files_s){ 
  load(i)
}

survival_windows <- bind_rows(
  makeWindowTable(ARTR_survival_monthly_ClimWin, 'ARTR', 'survival'), 
  makeWindowTable(HECO_survival_monthly_ClimWin, 'HECO', 'survival'), 
  makeWindowTable(POSE_survival_monthly_ClimWin, 'POSE', 'survival'), 
  makeWindowTable(PSSP_survival_monthly_ClimWin, 'PSSP', 'survival')
)
write_csv( survival_windows, 'output/survival_models/top_survival_windows.csv')

spList <- c('ARTR', 'HECO', 'POSE', 'PSSP')
sp <- 'ARTR'

for( sp in spList){ 
  make_window_plots(growth_windows, sp )
  make_window_plots(survival_windows, sp )
}


# Small plant growth 
result_files <- dir('output/growth_models/', pattern = '_small_plant_growth_mer_monthly_ClimWin.rda', full.names = T)

for( i in result_files){ 
  load(i)
}

growth_windows <- bind_rows( 
  makeWindowTable(ARTR__small_plant_growth_monthly_ClimWin, 'ARTR', 'growth'), 
  makeWindowTable(HECO__small_plant_growth_monthly_ClimWin, 'HECO', 'growth'), 
  makeWindowTable(POSE__small_plant_growth_monthly_ClimWin, 'POSE', 'growth'), 
  makeWindowTable(PSSP__small_plant_growth_monthly_ClimWin, 'PSSP', 'growth')

) 

write_csv( growth_windows, 'output/growth_models/top_small_growth_windows.csv')

spList <- c('ARTR', 'HECO', 'POSE', 'PSSP')
for( sp in spList){ 
  plot_small_plant_windows(growth_windows, sp )
}
