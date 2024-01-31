####################################################################################
#
# run scripts to fetch all demographic data 
# 
####################################################################################

rm(list = ls())

# growth
source("code/prepare_data/growth/fetchGrowthData.R")
source("code/prepare_data/growth/get_ARTR_growth.R")

# survival
source("code/prepare_data/survival/fetchSurvData.R")

get_growth_files <- dir('code', pattern = 'get.*growth\\.R', recursive = TRUE, full.names = TRUE)

get_survival_files <- dir('code', pattern = 'get.*survival\\.R', recursive = TRUE, full.names = TRUE)

# recruitment
get_recruitment_files <- dir('code', pattern = 'get.*recruitment\\.R', recursive = TRUE, full.names = TRUE)

lapply( c(get_survival_files, get_recruitment_files), source) 


#source('code/prepare_data/growth/clean_size_data.R')
