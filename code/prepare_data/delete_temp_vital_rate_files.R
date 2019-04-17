rm (list = ls()) 


combos <- expand.grid(spp = c('ARTR', 'HECO', 'POSE', 'PSSP'), vr = c('growth','survival', 'recruitment')) 

temp_files <- file.path( 'data/temp_data', paste0( combos$spp, '_', combos$vr, '.RDS'))

to_remove <- temp_files[ file.exists(temp_files) ] 

file.remove(to_remove)  # Clean up temporary files 
