rm(list = ls())
library(rstan)
library(tidyverse)

source('analysis/stan_data_functions.R')

top_model <- readRDS('output/stan_fits/ARTR_growth_model_3_top_model.RDS')
null_model <- readRDS('output/stan_fits/ARTR_growth_model_NULL_MOD_top_model.RDS')

test <- rstan::summary(top_model, 'hold_mu')
test$summary[, 1]

sum( get_lpd(top_model) )
sum( get_lpd(null_model))


top_model <- readRDS('output/stan_fits/POSE_growth_model_4_top_model.RDS')
null_model <- readRDS('output/stan_fits/POSE_growth_model_NULL_MOD_top_model.RDS')

sum(get_lpd(top_model))
sum(get_lpd(null_model))

top_model <- readRDS('output/stan_fits/PSSP_growth_model_4_top_model.RDS')
null_model <- readRDS('output/stan_fits/PSSP_growth_model_NULL_MOD_top_model.RDS')

sum(get_lpd(top_model))
sum(get_lpd(null_model))

top_model <- readRDS('output/stan_fits/ARTR_survival_model_6_top_model.RDS')
null_model <- readRDS('output/stan_fits/ARTR_survival_model_NULL_MOD_top_model.RDS')

sum( get_lpd(top_model) )
sum( get_lpd(null_model))



top_model <- readRDS('output/stan_fits/POSE_survival_model_3_top_model.RDS')
null_model <- readRDS('output/stan_fits/POSE_survival_model_NULL_MOD_top_model.RDS')

sum(get_lpd(top_model))
sum(get_lpd(null_model))

top_model <- readRDS('output/stan_fits/PSSP_survival_model_1_top_model.RDS')
null_model <- readRDS('output/stan_fits/PSSP_survival_model_NULL_MOD_top_model.RDS')

sum(get_lpd(top_model))
sum(get_lpd(null_model))



