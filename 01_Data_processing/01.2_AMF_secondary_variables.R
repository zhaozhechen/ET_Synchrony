# Author: Zhaozhe Chen
# Update date: 2025.8.11

# This code matches soil properties to AMF sites
# This code is done, no need to rerun, unless need to match a few sites with missing values

# ---------- Global ---------------
library(here)

# Input path to AMF site info
site_info <- read.csv(here("00_data","ameriflux_site_info.csv"))
# Input path to Soil properties
# Note: soil texture classification and code can be found in https://ftp.ems.psu.edu/pub/data/1995-0795/textclass.ascii
Soil_path <- here("00_Data")
# Soil hydraulic traits lookup table
Soil_lookup <- read.csv(here("00_Data","Soil_texture_lookup.csv"))

# Source functions for AMF processing
source(here("01_Data_processing","AMF_processing_functions.R")) 
  
# ------ Main ----------
# Match soil texture to AMF sites ======
# Note: only took the first layer
texture <- extract_nc(Soil_path,"Texture")

# Loop over all sites and get the soil texture for each site
# Initialize a vector to store soil texture for all sites
texture_all <- c()
for(i in 1:nrow(site_info)){
  # Get the coordinates at this site
  site_lon <- site_info$longitude[i]
  site_lat <- site_info$latitude[i]
  # Extract the nearest soil texture for this site
  site_texture <- get_nearest_value(site_lat = site_lat,
                                    site_lon = site_lon,
                                    grid_value = texture$layer,
                                    grid_lat = texture$lat,
                                    grid_lon = texture$lon,
                                    max_r = 1)
  texture_all <- c(texture_all,site_texture)
}
site_info$soil_texture <- texture_all
# Match corresponding soil hydraulic traits with the sites
site_info_update <- site_info %>%
  left_join(Soil_lookup,by=c("soil_texture" = "Class.code"))
write.csv(site_info_update,here("00_Data","ameriflux_site_info_update.csv"))
