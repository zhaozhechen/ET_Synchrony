# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.9.22

# This code is to extract Synchrony metrics (peak TE, memory, lag, aggTE)
# Peak daily TE: maximum significant TE normalized by Shannon entropy of the sink within the first 24 hours
# Lag: best lag corresponding to peak TE
# Daily aggregate TE: The total magnitude/contribution of predictive influence distributed across 24 hours, defined as the sum of TE for the set of lags with significant TE within 24 hours.
# Memory: Output 5 memory options:
# Option 1: From lag 0 to first non-significant lag
# Option 2: Width of continuous significant TE that includes peak TE
# Option 3: From peak to first non-significant lag after the peak
# Option 4: Total significant duration (cumulative hours)
# Option 5: From lag 0 to first lag after the peak

# Note: Focus on full TS, GS, and NGS for now. No years results included

# ------ Global ------ 
library(dplyr)
library(future)
library(future.apply)
library(progressr)

# Input path for TE_df
TE_df_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/TE_df/"
# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")
# Source functions
source("03_PN_construction/Synchrony_metrics_functions_v2.R")
source("05_Visualization/Plotting_functions.R")

max_lag <- 72 # Maximum lag to consider
# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")

# ------- Main ------
# All variable combinations
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# Sites to process
Site_IDs <- Site_info$site_id

syc_metrics_all_sites_df <- c()
# Process for all sites
for(i in 1:length(Site_IDs)){
  Site_ID <- Site_IDs[i]
  syc_metrics_all_sites <- syc_site(Site_ID,var_comb)
  syc_metrics_all_sites_df <- rbind(syc_metrics_all_sites_df,syc_metrics_all_sites)
  message(i)
}

# Output this combined df
write.csv(syc_metrics_all_sites_df,"03_PN_construction/Results/Syc_metrics_all_sites_5mem.csv")


