# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.8.17

# This code is to extract Synchrony metrics (peak TE, memory, lag)


# ------ Global ------ 
library(dplyr)

# Input path for TE_df
TE_df_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/TE_df/"
# Input path for Var_plots
Var_plots_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Var_plots/"
# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# ------- Main ------
# Get a list of Site_ID with valid TE output
Site_ls <- unique(substr(list.files(Var_plots_path),11,16))

setdiff(Site_info$site_id,Site_ls)

