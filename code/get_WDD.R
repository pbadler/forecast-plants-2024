# calculate wet degree days
library(dplyr)
library(lubridate)

VWC_threshold <- 10
start_month <- 2   # Sept is month 1 of the water year, March is 7

dayD <- read_csv('data/temp/daily_weather_for_models.csv')

# define function to find day at which 50% of cumulative wet degree days was reached
get_day50 <- function(x,y){
  target <- max(y)/2
  days <- which(y>target)
  day50 <- x[min(days)]
  return(day50)
}

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

# look at a time series
tmp <- subset(cumWDD,Treatment=="Control")
plot(tmp$date[1:2000],tmp$cumWDD[1:2000],type="l")

# aggregate to year
annual_WDD <- cumWDD %>% group_by(wat_year,Treatment) %>%
              summarize(WDD = max(cumWDD),
              day50 = get_day50(x=wat_day,y=cumWDD))

# test get_day50
x <- cumWDD$wat_day[cumWDD$wat_year==1950 & cumWDD$Treatment=="Control"]
y <- cumWDD$cumWDD[cumWDD$wat_year==1950 & cumWDD$Treatment=="Control"]
get_day50(x,y)
