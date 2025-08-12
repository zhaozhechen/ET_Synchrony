# Author: Zhaozhe Chen
# Date: 2025.8.11

# This code is to smooth LAI using Savitzky-Golay filter
# And calculate start of season (SOS) and end of season (EOS) for each year, at each sites


# --------- Global ----------
library(signal) # the function sgolayfilt is used to apply Savitzky-Golay smoothing filter

# Input path for LAI data
LAI_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/AMF_LAI/MCD15A3H_Lai_500m_"
# Input path to AMF site info
site_info <- read.csv(here("00_data","ameriflux_site_info_update.csv")) %>% select(-X)
# parameters to use in the SG filter
# window size
w <- 13
# degree of polynomial
d <- 4

# --------- Main ---------
for(i in 1:nrow(site_info)){
  Site_ID <- site_info$site_id[i]
  # Get LAI data
  LAI_df <- read.csv(paste0(LAI_path,Site_ID,".csv"))
  # Remove "F" and make it NA
  
  # Fill small gaps
  

  
  # (optional) fill small gaps so spacing stays regular
  lai_filled <- na.approx(lai, maxgap = 2, na.rm = FALSE)  # fill ≤2 consecutive gaps (~≤8 days)
  
  # Smooth
  lai_sg <- sgolayfilt(lai_filled, p = p, n = w, m = 0)
  # Smooth
  LAI_smooth <- sgolayfilt(LAI_df$Lai,p=d,n=w)
  
  
  # Calculate SOS and EOS
  
  #Make plots
  
  # Match with site_info
}


