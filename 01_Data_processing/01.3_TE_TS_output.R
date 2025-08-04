# Author: Zhaozhe Chen
# Update Date: 2025.8.4

# This code uses cleaned AMF hourly data as input,
# and outputs final time series of target variables for TE implementation for all sites
# One file for each site for simplicity

# ---------- Global ---------------
library(here)
library(dplyr)
library(future)
library(future.apply)
library(progressr)

# Input path to AMF site info, which also includes soil info
site_info <- read.csv(here("00_data","ameriflux_site_info_update.csv"))
# Input path for all cleaned hourly AMF dataset
AMF_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/02_AMF_cleaned/AMF_Hourly/"

# Source functions for AMF processing
source(here("01_Data_processing","AMF_processing_functions.R")) 

# Set parallel session
plan(multisession,workers = availableCores()-1)

# There are 5 cleaned AMF files, arrayid determines which file to process
arrayid <- 5

# ------ Main ----------
# Read in the cleaned AMF hourly data
AMF_df <- read.csv(paste0(AMF_path,"AMF_cleaned_hourly_",arrayid,".csv"))
# Calculate ET. convert unit from Wm-2 to mmday-1
AMF_df$ET <- AMF_df$LE_F * 3600 *24/(2.45*10^6)
AMF_df$ET[AMF_df$ET < 0] <- 0

# Get the Site ID in this file
Site_ID_ls <- unique(AMF_df$Site_ID)
# Loop over each site and calculate its target variables
progressr::with_progress({
  p <- progressor(along = 1:length(Site_ID_ls))
  
  future_lapply(1:length(Site_ID_ls),function(i){
    Site_ID <- Site_ID_ls[i]
    # Only keep data for this site
    Site_df <- AMF_df %>%
      filter(Site_ID == Site_ID)
    # Calculate soil water potential based on soil moisture and soil texture at this site
    # Reference: Lowman et al., 2023
    Site_df$psi_soil <- Cal_psisoil(Site_df,site_info,Site_ID)
    # Output this df
    write.csv(Site_df,paste0(AMF_path,"AMF_sites_hourly/AMF_hourly_",Site_ID,".csv"))
    p()
  })
})








