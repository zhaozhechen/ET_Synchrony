# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.8.17

# This code is to extract Synchrony metrics (peak TE, memory, lag)
# Peak TE: maximum significant TE normalized by Shannon entropy of the sink
# Memory: time required for TE to become insignificant
# Lag: best lag corresponding to peak TE
# Note: Focus on full TE, GS, and NGS for now. No years results included

# ------ Global ------ 
library(dplyr)

# Input path for TE_df
TE_df_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/TE_df/"
# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")
# Source functions
source("03_PN_construction/Synchrony_metrics_functions.R")

max_lag <- 72 # Maximum lag to consider
# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")

# ------- Main ------
# Determines which Site to proceed
arrayid <- 1
Site_ID <- Site_info$site_id[arrayid]
# All variable combinations
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# For full year =======
file_name <- paste0("TE_df_ls_full_TS_",Site_ID,".rds")


















