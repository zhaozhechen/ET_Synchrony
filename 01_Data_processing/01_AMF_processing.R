# Author: Zhaozhe Chen
# Update date: 2025.6.30

# This code conducts QC for AMF dataset
# Output SM, ET
# More descriptions here (need update)
# This is a test code, just to get hourly SM and ET for US-Ne1




# --------- Global ---------
library(here)
library(dplyr)
library(lubridate)

# Input path to AMF dataset
Data_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/"
# Source functions for AMF processing
source(here("01_Data_processing","AMF_processing_functions.R"))
# Output path for figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_AMF_TS/"
# The raw hourly dataset is divided into 5 files for all sites, so use
# arrayid to determine which file to process (1-5)
arrayid <- 2
# This is a list of required variables for output
var_ls <- c("SWC", # Soil Water Content, Unit: %
            "TA", # Air temperature, Unit: degree C
            "VPD", # Vapor pressure deficit, Unit: hPa
            "P_F", # Precipitation (to be separated with pressure), Unit: mm
            "WS", # Wind speed, Unit: ms-1
            "USTAR", # Friction velocity, Unit: ms-1
            "NETRAD", # Net radiation, Unit: Wm-2
            "TS", # Soil temperature, Unit: degree C
            "LE_F" # Latent heat flux, Unit: Wm-2
            )

# ---------- Main ---------
gc()
# Get site info of each site
#Site_info <- read.csv(paste0(Data_path,"01_AMF_raw/ameriflux_site_info.csv"))

# Get file names for the 5 files
file_names <- list.files(paste0(Data_path,"01_AMF_raw"),pattern = "_data_",full.names = TRUE)
# Read in raw AMF dataset
AMF_raw <- read.csv(file_names[arrayid])
# For each required variable, QC all sub-variables for them, and get the average of all sub-variables for each variable
AMF_df <- Var_QC_all(var_ls,AMF_raw)

# Aggregate all data to hourly data (some AMF sites have half-hourly data)
AMF_df <- AMF_df %>%
  mutate(Time_hour = floor_date(Time,unit="hour")) %>%
  group_by(Time_hour,Site_ID) %>%
  summarise(across(where(is.numeric),~ mean(.x,na.rm=TRUE)),
            .groups = "drop") %>%
  rename(Time = Time_hour)
  



# Only keep data for the target site
AMF_site <- AMF_raw %>%
  filter(site_id == "US-Ne1")
rm(AMF_raw)
# Extract SWC data
SWC <- Mean_QC_Var("SWC",AMF_site,Coarse = FALSE)
# Get SM Unit: m3/m3
SM <- SWC/100
# Get ET Unit: Convert from W/m2 to mm/day
ET <- Mean_QC_Var("LE_F",AMF_site,Coarse = FALSE) * 3600 * 24/(2.45*10^6)
ET[ET < 0] <- 0


# Output this data frame
write.csv(df_hourly,paste0(Output_path,"US-Ne1_hourly_SM_ET_test.csv"))




