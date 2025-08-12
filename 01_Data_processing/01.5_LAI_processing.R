# Author: Zhaozhe Chen
# Date: 2025.8.12

# This code is to smooth LAI using Savitzky-Golay filter
# And calculate start of season (SOS) and end of season (EOS) for each year, at each sites

# --------- Global ----------
library(signal) # the function sgolayfilt is used to apply Savitzky-Golay smoothing filter
library(dplyr)
library(here)
library(zoo) # the function na.approx is used to replace NA by interpolation

# Input path for LAI data
LAI_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/AMF_LAI/MCD15A3H_Lai_500m_"
# Input path to AMF site info
site_info <- read.csv(here("00_data","ameriflux_site_info_update.csv")) %>% select(-X)
# parameters to use in the SG filter
# window size
windowsize <- 13
# degree of polynomial
degree_p <- 4

# --------- Main ---------
for(i in 1:nrow(site_info)){
  Site_ID <- site_info$site_id[i]
  # Get LAI data
  LAI_df <- read.csv(paste0(LAI_path,Site_ID,".csv"))
  LAI_df <- LAI_df %>%
    mutate(
      # Remove "F" and make it NA
      Lai = as.numeric(if_else(Lai == "F",NA,Lai)),
      # Fill small gaps, fill <=2 consecutive gaps with NA
      LAI_filled = na.approx(Lai,maxgap=2,na.rm=FALSE),
      # Apply Savistzky-Golay smoothing filter
      LAI_smoothed = sgolayfilt(LAI_filled,p = degree_p,n = windowsize)
    )
  # Calculate LAI thresholds (50% of max and min LAI)
  
  # Calculate SOS and EOS
  
  #Make plots
  
  # Match with site_info
}


