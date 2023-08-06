rm(list = ls())
library(tidyverse)
library(lubridate)
library(climwin)
library(lme4)
library(optimx)

# LMER optimization options
control_lmer = lmerControl(
  optimizer = "optimx",
  calc.derivs = FALSE,
  optCtrl = list(
    method = "nlminb",
    starttests = FALSE,
    kkt = FALSE
  )
)
control_lmer$optCtrl$eval.max <- 1e8
control_lmer$optCtrl$iter.max <- 1e8

source('code/analysis/functions.R')

# Variables -------------------------------- : 
last_year <- 2010 # last year of training data, everything before and including this year is used 
sp_list <- c('ARTR') #, 'HECO', 'POSE', 'PSSP')

# ClimWin Window Settings Monthly
window_open_max <- 24
window_open_min <- 1
window_exclude_dur <- 1
window_exclude_max <- 5

size_cutoff <- -1 # size division for big/small 

load('output/tempDataForTesting.rda')

newVars <- as.list( daily_weather[, addVars_list$addVars])

growthWin2 <- slidingwin(xvar = newVars,
                         cdate = daily_weather$date_reformat,
                         bdate = addVars_list$data2$date_reformat,
                         baseline = m_baseline,
                         cinterval = 'month',
                         range = c(window_open_max, window_open_min),
                         #exclude = c(window_exclude_dur, window_exclude_max),
                         type = "absolute", refday = c(15, 06),
                         stat = 'mean',
                         func = c('lin'),
                         cv_by_cohort = TRUE, ncores = 8)

load('output/tempDataForTesting2.rda')
baseline
modeloutput


formula( baseline)
cross_validate(fold = 1, modeldat = modeldat, baseline = baseline, modeloutput = modeloutput)
