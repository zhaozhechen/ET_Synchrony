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

# Input path to cleaned AMF dataset
AMF_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/02_AMF_cleaned/AMF_Hourly/"
# Input path to Soil properties
Soil_path <- here("00_Data")
# Source functions for AMF processing
source(here("01_Data_processing","AMF_processing_functions.R")) 
  
# ------ Main ----------
arrayid <- 1
# Get file names for the hourly AMF dataset
file_names <- list.files(AMF_path,full.names = TRUE)
# Read in hourly AMF dataset
AMF_df <- read.csv(file_names[arrayid])
# Get Soil property file paths
# Soil porosity
Soil_Property_file <- list.files(Soil_path,full.names = TRUE)[grepl("Porosity",list.files(Soil_path))]
# 


#Soil_Property <- stack(paste0(Soil_path,"CON"))






