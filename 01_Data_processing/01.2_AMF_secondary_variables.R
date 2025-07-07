# Author: Zhaozhe Chen
# Update date: 2025.7.6

# This code is to calculate secondary variables from cleaned AmeriFlux dataset
# The calculated secondary variables include:
# Soil water potential
# PET
# 
# ---------- Global ---------------
library(dplyr)
library(here)
library(ncdf4)

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
test <- extract_nc(Soil_path,"Texture")






# Get file names for the hourly AMF dataset
#file_names <- list.files(AMF_path,full.names = TRUE)
# Read in hourly AMF dataset
#AMF_df <- read.csv(file_names[arrayid])