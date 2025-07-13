# Author: Zhaozhe Chen
# Update Date: 2025.7.7

# This code is to output final time series of target variables for TE implementation
# Next steps:
# Calculate PET,etc.


# ---------- Global ---------------
library(here)

# Input path to AMF site info
site_info <- read.csv(here("00_data","ameriflux_site_info_update.csv"))
# Source functions for AMF processing
source(here("01_Data_processing","AMF_processing_functions.R")) 
# Hourly SM and ET data at US-Ne1. This is test dataset
AMF_df <- read.csv(here("00_Data","AMF_hourly_test.csv"))
Site_ID <- "US-Ne1"

# ------ Main ----------
# Calculate soil water potential based on soil moisture and soil texture
AMF_df$psi_soil <- Cal_psisoil(AMF_df,site_info,Site_ID)

# Calculate ET convert unit from Wm-2 to mmday-1
AMF_df$ET <- AMF_df$LE_F * 3600 *24/(2.45*10^6)
AMF_df$ET[AMF_df$ET < 0] <- 0

# Output this df
write.csv(AMF_df,here("00_Data","Data_processed",paste0("AMF_hourly_",Site_ID,".csv")))






