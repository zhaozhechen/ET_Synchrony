# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Update date: 2025.6.24

# ------- Global --------
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(RTransferEntropy)
library(entropy)
library(future)
library(lubridate)
library(progress)

source(here("02_TE_implementation","TE_implementation_functions.R"))
source(here("05_Visualization","Plotting_functions.R"))

# Cleaned daily AMF dataset
AMF_df <- read.csv(here("00_Data","AMF_daily_all_sites.csv"))
# Output path for figures
Output_path <- here("02_TE_implementation","Results/")

# Number of bins for TE discritization
n_bin <- 11
my_seed <- 111
max_lag <- 75
# Set parallel session
plan(multisession,workers = availableCores()-1)

# Site ID to test
Site_test <- "US-Ne1"

# -------- Main ---------
# Get data for this AMF site
AMF_df <- AMF_df %>%
  filter(Site_ID == Site_test)

# Get change in SM, VPD, and ET
delta_SM <- delta_TS(AMF_df,"SM")
delta_VPD <- delta_TS(AMF_df,"VPD")
delta_ET <- delta_TS(AMF_df,"ET")

# Visualize the time series
# Time for the delta TS
TS_time <- as.Date(tail(AMF_df$Date,-1))
# Output plots for all these delta TS
print_g(TS_plot(delta_SM,TS_time,bquote(Delta~SM)),
        paste0("Delta_SM_TS_",Site_test),
        8,6)
print_g(TS_plot(delta_VPD,TS_time,bquote(Delta~VPD)),
        paste0("Delta_VPD_TS_",Site_test),
        8,6)
print_g(TS_plot(delta_ET,TS_time,bquote(Delta~ET)),
        paste0("Delta_ET_TS_",Site_test),
        8,6)

# Run TE from delta_SM to delta_ET
SM_TE_df <- Cal_TE(delta_SM,delta_ET,max_lag = max_lag)
VPD_TE_df <- Cal_TE(delta_VPD,delta_ET,max_lag = max_lag)

# Plot TE vs lag
g_SM_TE <- TE_lag_plot(SM_TE_df,bquote(Delta~SM~"->"~Delta~ET))
g_VPD_TE <- TE_lag_plot(VPD_TE_df,bquote(Delta~VPD~"->"~Delta~ET))
# Combine these two plots
g_TE <- plot_grid(g_SM_TE,g_VPD_TE,nrow=1,align="hv")
print_g(g_TE,paste0("TE_",Site_test),6,3)

