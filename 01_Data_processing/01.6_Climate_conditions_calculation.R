# Author: Zhaozhe Chen
# Update Date: 2025.9.22

# ---------- Global ---------------
library(here)
library(dplyr)


# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")
# Input path for all cleaned hourly AMF dataset
AMF_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/02_AMF_cleaned/AMF_Hourly/AMF_sites_hourly_update/"

# Source functions for AMF processing
source(here("01_Data_processing","AMF_processing_functions.R")) 

# Manually add CVM, which has the same gs as CRO
# Manually add BSV
PFTlist <- c('ENF','EBF','DNF','DBF','MF','CSH','OSH','WSA','SAV','GRA','CRO','CVM','BSV') # plant functional type
GSlist <- 1/c(125,150,150,100,125,300,170,300,70,40,40,40,1000)# m/s, corresponding stomatal conductance from  MODIFIED_IGBP_MODIS_NOAH
# Parameters for PET calculation Unit mm/day
rhoa  <- 1.225  # kg/m3
Cp    <- 1005   # J/kg/K
gamma <- 66     # Pa/K
Lv    <- 2453e6 # J/m3

# ------------- Main ----------

# Initialize vector to store output
AI_all <- c()
ER_all <- c()

for(i in 1:nrow(Site_info)){
  Site_ID <- Site_info$site_id[i]
  Site_PFT <- Site_info$IGBP_veg[i]
  # Get gs
  gs <- GSlist[PFTlist==Site_PFT]  
  # Read in AMF_df
  AMF_df <- read.csv(paste0(AMF_path,"AMF_hourly_",Site_ID,".csv")) %>%
    select(-c(1,2))
  
  AMF_df <- AMF_df %>%
    # Calculate Delta Unit Pa/K
    mutate(DEL = PM_delta(TA + 273.15),
           ga = 1/(WS/USTAR^2 + 6.2*USTAR^(-2/3))) %>%
    # Get hourly PET mm/hour
    mutate(hourly_PET = (DEL*NETRAD + rhoa*Cp*ga*VPD*1000)/(DEL+gamma*(1+ga/gs))/Lv*1000*3600)
  
  # Aridity index calculation --------
  # Keep only rows when P and PET both exisit
  AI_rows <- AMF_df %>%
    filter(!is.na(P_F),!is.na(hourly_PET))
  # Get total P and total PET
  total_P <- sum(AI_rows$P_F)
  total_PET <- sum(AI_rows$hourly_PET)
  AI <- total_P/total_PET

  # Evaporative ratio (ER) calculation -------
  # Keep only rows when P and ET both exisit
  ER_rows <- AMF_df %>%
    filter(!is.na(P_F),!is.na(ET))
  total_P <- sum(ER_rows$P_F)
  total_ET <- sum(ER_rows$ET)
  ER <- total_P/total_ET
  
  AI_all <- c(AI_all,AI)
  ER_all <- c(ER_all,ER)
  
  message(i)

}

# Put them into df
AI_ER_df <- data.frame(Site_ID = Site_info$site_id,
                       AI = AI_all,
                       ER = ER_all)

# Output this df
write.csv(AI_ER_df,"00_Data/ameriflux_site_AI_ER.csv")



