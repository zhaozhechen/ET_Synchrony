# Author: Zhaozhe Chen
# Update date: 2025.7.7

# This code is to run TE at hourly scale

# -------- Global ------------
library(here)
# Source functions for data processing, TE implementation, and visualization
source(here("01_Data_processing","AMF_processing_functions.R"))
source(here("02_TE_implementation","TE_implementation_functions.R"))
source(here("05_Visualization","Plotting_functions.R"))

# Input path for hourly AMF data
# Note: this is a test dataset at US-Ne1, need to be updated later
AMF_df <- read.csv(here("00_Data","Data_processed","AMF_hourly_US-Ne1.csv"))
# Output path for the results
Output_path <- here("02_TE_implementation","Results","Hourly_TE_test_US-Ne1")
# Test at US-Ne1
Site_ID <- "US-Ne1"

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

# These are folding parameters to deal with extreme values (outliers) in the time series
# i.e., extreme values will be binned into the first or last bin
lower_qt <- 0.001
upper_qt <- 0.999

# ----------- Main -----------
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

# Put them in a df
df <- data.frame(
  Time = AMF_df$Time[2:nrow(AMF_df)],
  delta_log10_psi_soil,
  delta_VPD,
  delta_ET)

# Get diurnal anomaly
df$delta_log10_psi_soil_anomaly <- 
  Cal_diurnal_anomaly(df,"delta_log10_psi_soil",5)
df$delta_ET_anomaly <-
  Cal_diurnal_anomaly(df,"delta_ET",5)
df$delta_VPD_anomaly <- 
  Cal_diurnal_anomaly(df,"delta_VPD",5)

# TE from delta_VPD_anomaly -> delta_ET_anomaly
# Timing the TE calculation
start_time <- Sys.time()
VPD_TE_df <- Cal_TE_main(var1 = df$delta_VPD_anomaly,
                     var2 = df$delta_ET_anomaly,
                     max_lag = max_lag,
                     nbins = n_bin,
                     alpha = alpha,
                     nshuffle = nshuffle,
                     upper_qt = upper_qt,
                     lower_qt = lower_qt,
                     ZFlag_Source = FALSE,
                     ZFlag_Sink = FALSE)
end_time <- Sys.time()
run_time <- as.character(end_time - start_time)
writeLines(run_time,con=paste0(Output_path,"/TE_VPD_log.txt"))
write.csv(VPD_TE_df,paste0(Output_path,"/TE_df_US-Ne1_hourly_VPD_ET.csv"))
g <- TE_lag_plot(VPD_TE_df,"VPD->ET","None")
print_g(g,"TE_US-Ne1_hourly_VPD_ET",6,4)

# TE from delta_log10_psi_soil_anomaly -> delta_ET_anomaly
# Timing the TE calculation
start_time <- Sys.time()
psi_TE_df <- Cal_TE_main(var1 = df$delta_log10_psi_soil_anomaly,
                     var2 = df$delta_ET_anomaly,
                     max_lag = max_lag,
                     nbins = n_bin,
                     alpha = alpha,
                     nshuffle = nshuffle,
                     upper_qt = upper_qt,
                     lower_qt = lower_qt,
                     ZFlag_Source = FALSE,
                     ZFlag_Sink = FALSE)
end_time <- Sys.time()
run_time <- as.character(end_time - start_time)
writeLines(run_time,con=paste0(Output_path,"/TE_psi_log.txt"))
write.csv(psi_TE_df,paste0(Output_path,"/TE_df_US-Ne1_hourly_psi_ET.csv"))
g <- TE_lag_plot(psi_TE_df,"psi->ET","None")
print_g(g,"TE_US-Ne1_hourly_psi_ET",6,4)

# Plot TS of all variables
# For psi_soil
g1 <- TS_plot(AMF_df$psi_soil,AMF_df$Time,bquote(psi[soil]~"("~kPa~")"))
print_g(g1,paste0(Site_ID,"_psi_soil_TS"),
        10,6)
# For log10(psi_soil)
g2 <- TS_plot(AMF_df$log10_psi_soil,AMF_df$Time,bquote(log10~psi[soil]~"("~kPa~")"))
print_g(g2,paste0(Site_ID,"_log10_psi_soil_TS"),
        10,6)
# For delta_log10_psi_soil
g3 <- TS_plot(df$delta_log10_psi_soil,df$Time,bquote(Delta~log10~psi[soil]~"("~kPa~")"))
print_g(g3,paste0(Site_ID,"_delta_log10_psi_soil_TS"),
        10,6)
# For ET
g4 <- TS_plot(AMF_df$ET,AMF_df$Time,"ET (mm/day)")
print_g(g4,paste0(Site_ID,"_ET_TS"),
        10,6)
# For delta_ET
g5 <- TS_plot(df$delta_ET,df$Time,bquote(Delta~ET~"("~mm/day~")"))
print_g(g5,paste0(Site_ID,"_delta_ET_TS"),
        10,6)
# For delta_log10_psi_soil_anomaly
g6 <- TS_plot(df$delta_log10_psi_soil_anomaly,df$Time,bquote(Delta~log10~psi[soil]~anomaly~"("~kPa~")"))
print_g(g6,paste0(Site_ID,"_delta_log10_psi_soil_anomaly_TS"),
        10,6)
# For delta_ET_anomaly
g7 <- TS_plot(df$delta_ET_anomaly,df$Time,bquote(Delta~ET~anomaly~"("~mm/day~")"))
print_g(g7,paste0(Site_ID,"_delta_ET_anomaly_TS"),
        10,6)
# For VPD
g8 <- TS_plot(AMF_df$VPD,AMF_df$Time,"VPD (kPa)")
print_g(g8,paste0(Site_ID,"_VPD_TS"),
        10,6)
# For delta_VPD
g9 <- TS_plot(df$delta_VPD,df$Time,bquote(Delta~VPD~"("~kPa~")"))
print_g(g9,paste0(Site_ID,"_delta_VPD_TS"),
        10,6)
# For delta_VPD_anomaly
g10 <- TS_plot(df$delta_VPD_anomaly,df$Time,bquote(Delta~VPD~anomaly~"("~kPa~")"))
print_g(g10,paste0(Site_ID,"_delta_VPD_anomaly_TS"),
        10,6)
