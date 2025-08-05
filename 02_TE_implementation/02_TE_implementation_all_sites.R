# Author: Zhaozhe Chen
# Update date: 2025.8.4

# This code is to run TE at hourly scale for all AMF sites
# Calculate TE between variable pairs of:
# Soil water potential, ET, and T
# For both directions

# Output includes:
# 

# -------- Global ------------
library(here)
# Source functions for data processing, TE implementation, and visualization
source(here("01_Data_processing","AMF_processing_functions.R"))
source(here("02_TE_implementation","TE_core_codes.R"))
source(here("05_Visualization","Plotting_functions.R"))

# Input path for hourly AMF data
AMF_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/02_AMF_cleaned/AMF_Hourly/AMF_sites_hourly/"
# Input path to AMF site info, which also includes soil info
site_info <- read.csv(here("00_data","ameriflux_site_info_update.csv"))
# Output path for the results
Output_path <- here("02_TE_implementation","Results","Hourly_TE_all_sites")

# Parameters for TE implementation
n_bin <- 11 # Number of bins for TE discritization of continuous data (e.g., SM)
max_lag <- 72 # Maximum lag to consider (This should be adjusted according to the processes and the temporal resolution of data)
Lag_Dependent_Crit <- TRUE # Determine if critical TE is lag-dependent
nshuffle <- 300 # Number of shuffles (bootstrap) for critical TE for statistical inference
alpha <- 0.05 # Confidence level for critical TE
# Set parallel session
plan(multisession,workers = availableCores()-1)
# Ensure reproducibility
set.seed(111)
# Determines if zero should be adjusted for the Sink and Source variables
ZFlagSink <- FALSE
ZFlagSource <- FALSE

# Whether to plot TS and histogram of the source and sink
TS_Hist_plot <- TRUE

# These are folding parameters to deal with extreme values (outliers) in the time series
# i.e., extreme values will be binned into the first or last bin
lower_qt <- 0.001
upper_qt <- 1-lower_qt

# 2 colors for growing season and non-growing season
season_color <- brewer.pal(3,"Set2")[1:2]

# Determines which site to process
arrayid <- 73

# -------- Main ----------
Site_ID <- site_info$site_id[arrayid]

# Step 1. Data processing and Calculation of secondary variables ---------------- 
# Read in hourly data for this site
AMF_df <- read.csv(paste0(AMF_path,"AMF_hourly_",Site_ID,".csv"))
# Standardize the time of the df
AMF_df <- Standardize_time(AMF_df)
# Convert psi_soil to log(psi_soil) to reduce skewness
AMF_df$log10_psi_soil <- log10(AMF_df$psi_soil)
# Calculate change in delta_log10_psi_soil
delta_log10_psi_soil <- delta_TS(AMF_df,"log10_psi_soil")
# Calculate delta_ET
delta_ET <- delta_TS(AMF_df,"ET")
# Calculate delta_VPD
delta_VPD <- delta_TS(AMF_df,"VPD")
# Calculate delta_TA
delta_TA <- delta_TS(AMF_df,"TA")

# Put them in a df
df <- data.frame(
  Time = AMF_df$Time[2:nrow(AMF_df)],
  delta_log10_psi_soil,
  delta_VPD,
  delta_ET,
  delta_TA)

# Get diurnal anomaly
df<- Cal_diurnal_anomaly(df,"delta_log10_psi_soil",5)
df<- Cal_diurnal_anomaly(df,"delta_VPD",5)
df<- Cal_diurnal_anomaly(df,"delta_ET",5)
df <- Cal_diurnal_anomaly(df,"delta_TA",5)

# Also add precipitation as the original data
df$PPT <- AMF_df$P_F[2:nrow(AMF_df)]

# Define growing season (GS)
# Define GS as May to Sep, could be revised if needed
df <- df %>%
  mutate(
    GS =  if_else(month(Time) %in% 5:9,"GS","Non-GS")
  )

# Step 2. Make TS and histogram plots of input variables ----------------------
# Get TS plots and and distribution for ET, psi, VPD, TA, and P
# Note: all the first four variables are diurnal anomaly of delta_TS, only P is the original daily P
g_psi <- var_plot_TS_Hist("delta_log10_psi_soil_anomaly",y_title = bquote(Delta~psi),
                          df,my_color = season_color,ZFlag = FALSE,nbins=n_bin)

# Add Site ID to the top of the 



