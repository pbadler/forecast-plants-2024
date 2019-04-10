#####################################################################################
#
# Make climate variables from seasonal climate 
#
#####################################################################################

rm(list = ls() )

# ------- load files ------------------------------------------------------------------ 

seasonal_clim <- readRDS('data/temp_data/seasonal_climate.RDS')
seasonal_VWC  <- readRDS('data/temp_data/seasonal_VWC.RDS')

# ------ calculate seasonal lags -----------------------------------------------------# 
#
#   Variable names follow these conventions: 
#   
#     First letter gives variable type: 
#       "P" is cumulative precipitation
#       "T" is average mean monthly temperature
#
#     Letters after the first period give the season aggregation window:
#
#       w  = winter (Q1) 
#       sp = spring (Q2)
#       su = summer (Q3)
#       f  = fall   (Q4)
#       a  = annual (Q1-4)
#
#     e.g. "P.sp" is the cumulative precipitation of the spring season and "P.w.sp" is
#     the cumulative precipitation of the winter and spring. 
#   
#     Number after the second period indicates the year of the transition, For example, 
#     "P.sp.0" gives the cumulative precipitation of year preceding the transition. 
#     Whereas "T.f.w.1" gives the average temperature of the fall and winter 
#     preceding the second year of the transition. "0" refers to year before first year, 
#     i.e. "lag effect" (sensu Adler). 
#
# -------------------------------------------------------------------------------------# 
seasonal_VWC$season <- factor( seasonal_VWC$season, c('winter', 'spring', 'summer', 'fall'), ordered = T)
seasonal_clim$season <- factor( seasonal_clim$season, c('winter', 'spring', 'summer', 'fall'), ordered = T)

q_VWC <- 
  seasonal_VWC %>% 
  spread(Treatment, avg) %>% 
  mutate( Drought = ifelse( year < 2012, Control, Drought  ),         # assign treatments prior to 2012 with Control level
          Irrigation = ifelse(year < 2012, Control, Irrigation)) %>% 
  gather( Treatment, avg, Control:Irrigation) %>% 
  ungroup() %>% 
  group_by(Treatment) %>% 
  arrange(Treatment, year, season) %>%
  mutate(VWC.sp.1 = avg, 
         VWC.sp.0 = lag( VWC.sp.1, 4),
         VWC.sp.l = lag( VWC.sp.0, 4),
         #VWC.sp.su.1 = rollapply(avg, 2, 'mean', na.rm = TRUE, align = 'right', fill = NA),
         #VWC.sp.su.0 = lag(VWC.sp.su.1, 4),
         #VWC.sp.su.l = lag(VWC.sp.su.0, 4),
         VWC.su.1 = lag(avg, 3),
         VWC.su.0 = lag(VWC.su.1, 4),
         VWC.su.l = lag(VWC.su.0, 4),
         VWC.f.1  = lag(avg, 2), 
         VWC.f.0  = lag(VWC.f.1, 4),
         VWC.f.l  = lag(VWC.f.0, 4)) %>%
         #VWC.a.1 = rollapply(avg, 4, 'mean', na.rm = TRUE, align = 'right', fill = NA), 
         #VWC.a.0 = lag( VWC.a.1,4), 
         #VWC.a.l = lag( VWC.a.0,4)) %>%
  filter( season == 'spring') %>% # plants are measured at the end of spring each year 
  dplyr::select( Treatment, Period, year, season, starts_with("VWC")) %>%
  ungroup() %>% 
  gather( var, val, starts_with('VWC')) %>% 
  filter( !is.na(val)) %>%
  spread( var, val) 

q_precip <- 
  seasonal_clim %>% 
  filter( var == 'PRCP_ttl') %>%
  group_by(Treatment) %>% 
  arrange(Treatment, year, season) %>%
  mutate(P.f.w.sp.1 = rollsum(val, 3, align = 'right', fill = NA), 
         P.f.w.sp.0 = lag(P.f.w.sp.1, 4),
         P.f.w.sp.l = lag(P.f.w.sp.0, 4),
         P.a.1      = rollsum(val, 4, align = 'right', fill = NA),
         P.a.0      = lag(P.a.1, 4),
         P.a.l  = lag(P.a.0, 4),
         P.su.1 = lag(val, 3),                 
         P.su.0 = lag(P.su.1, 4), 
         P.su.l = lag(P.su.0, 4)) %>% 
  filter( season == 'spring') %>% # plants are measured at the end of spring each year 
  dplyr::select( Treatment, Period, year, season, starts_with("P"))

q_temp <- 
  seasonal_clim %>% 
  filter( var == 'TAVG_avg' ) %>% 
  group_by(Treatment) %>% 
  arrange(Treatment, year, season) %>% 
  mutate( T.sp.1 = val, 
          T.sp.0 = lag(T.sp.1, 4),
          T.sp.l = lag(T.sp.0, 4), 
          T.su.1 = lag(val, 3), 
          T.su.0 = lag(T.su.1, 4), 
          T.su.l = lag(T.su.0, 4),
          T.f.1 = lag(val, 2), 
          T.f.0 = lag(T.f.1, 4), 
          T.f.l = lag(T.f.0, 4), 
          T.w.1 = lag(val, 1), 
          T.w.0 = lag(T.w.1, 4), 
          T.w.l = lag(T.w.0, 4)) %>% 
  filter( season == 'spring') %>% 
  dplyr::select( Treatment, Period, year, season, starts_with("T."))

allClim <- 
  q_precip %>% 
  left_join ( q_temp, by = c('Treatment', 'Period', 'season', 'year')) %>% 
  right_join ( q_VWC, by = c('Treatment', 'Period', 'season', 'year')) %>% 
  arrange( Treatment, year) 

NA_index <- apply( allClim, 1 , function(x) any(is.na(x)))

allClim <- allClim[-NA_index, ]

# adjust years so that they match the demographic data sets ------------------------------------------------------------------# 

allClim$year <- allClim$year - 1 # adjust to match assignment of year 0 as the reference year in demographic data sets

# ----------------------------------------------------------------------------------------------------------------------------# 

# -- calculate interactions --------------------------------------------------------------#

allClim <- data.frame(allClim)

# include interactions between temperature and VWC for each season
VWCvars <- which(is.element(substring(names(allClim),1,2),"VW")==T)
tmp_dat <- matrix(NA,nrow(allClim),length(VWCvars))
colnames(tmp_dat) <- paste0("VWCxT.",substring(names(allClim)[VWCvars],5))
for(i in 1:length(VWCvars)){
  do_season <- substring(names(allClim)[VWCvars[i]],5)
  iT <- which(names(allClim)==paste0("T.",do_season))
  tmp_dat[,i] <- allClim[,VWCvars[i]]*allClim[,iT]
}
allClim <- cbind(allClim,tmp_dat)

# ---- output ----------------------------------------------------------------------------# 

saveRDS( allClim , 'data/temp_data/all_clim_covs.RDS')
