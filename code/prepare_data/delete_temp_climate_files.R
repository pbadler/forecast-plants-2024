rm (list = ls()) 

rainfall_temp_file  <- 'data/temp/daily_station_dat_rainfall.RDS'
spot_temp_file <- 'data/temp/spotVWC.RDS'
swTreatment_temp_file <- 'data/temp/daily_swVWC_treatments.RDS'
decagon_data_temp_file <- 'data/temp/decagon_data_with_station_data.RDS'

file.remove(rainfall_temp_file,spot_temp_file, swTreatment_temp_file, decagon_data_temp_file )
