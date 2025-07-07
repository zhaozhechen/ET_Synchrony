# Author: Zhaozhe Chen
# Update date: 2025.7.6

# This code is to calculate secondary variables from cleaned AmeriFlux dataset
# The calculated secondary variables include:
# Soil water potential
# PET
# 
# ---------- Global ---------------
library(here)

# Input path to AMF site info
site_info <- read.csv(here("00_data","ameriflux_site_info.csv"))
# Input path to cleaned AMF dataset
AMF_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/02_AMF_cleaned/AMF_Hourly/"
# Input path to Soil properties
# Note: soil texture classification and code can be found in https://ftp.ems.psu.edu/pub/data/1995-0795/textclass.ascii
Soil_path <- here("00_Data")
# Source functions for AMF processing
source(here("01_Data_processing","AMF_processing_functions.R")) 
  
# ------ Main ----------
# Extract soil texture from nc files
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





# Get file names for the hourly AMF dataset
#file_names <- list.files(AMF_path,full.names = TRUE)
# Read in hourly AMF dataset
#AMF_df <- read.csv(file_names[arrayid])