# functions to calculate wet degree days

# define function to find day at which 50% of cumulative wet degree days was reached
get_day50 <- function(x,y){
  target <- max(y)/2
  days <- which(y>target)
  day50 <- x[min(days)]
  return(day50)
}

get_WDD <- function(weather, VWC_threshold, start_month){
  
  dayD <- weather #rename daily weather object
  
  # toss unnecessary columns
  dayD <- dplyr::select(dayD, date,Treatment,VWC,TAVG)
  
  # extract year, month and day
  dayD$year <- year(dayD$date)
  dayD$month <- month(dayD$date)
  dayD$jday <- yday(dayD$date)
  
  # convert to water year, month and date
  dayD$wat_year <- ifelse(dayD$month > 8, dayD$year + 1, dayD$year )
  dayD$wat_month <- ifelse(dayD$month > 8, dayD$month - 8, dayD$month + 4 )
  dayD$wat_day <- ifelse(dayD$jday > 243, dayD$jday - 243, dayD$jday + 122 )
  
  # calculate individual wet degree days
  dayD$is_wet <- dayD$VWC > VWC_threshold & dayD$TAVG > 0
  dayD$WDD <- ifelse(dayD$is_wet=="TRUE",dayD$TAVG,0)
  
  # calculate cumulative sum for each water year
  cumWDD <- dayD %>% filter(wat_month >= start_month & wat_month < 11) %>%
    group_by(wat_year,Treatment) %>% mutate(cumWDD = cumsum(WDD))
  
  # aggregate to year
  annual_WDD <- cumWDD %>% group_by(wat_year,Treatment) %>%
    summarize(WDD = max(cumWDD),
              day50 = get_day50(x=wat_day,y=cumWDD))
  
  return(annual_WDD)
  
}


