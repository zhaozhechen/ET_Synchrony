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
# All variable combinations
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# Loop over all sites
# Initialize a vector
syc_metrics_all_sites <- c()
for(arrayid in 1:nrow(Site_info)){
  Site_ID <- Site_info$site_id[arrayid]
  # For full year
  file_name_full_TS <- paste0("TE_df_ls_full_TS_",Site_ID,".rds")
  syc_metrics_full_TS <- cal_syc_metrics_all_pairs(file_name_full_TS)
  # For GS
  file_name_GS <- paste0("TE_df_ls_GS_",Site_ID,".rds")
  syc_metrics_GS <- cal_syc_metrics_all_pairs(file_name_GS)
  # For NGS
  file_name_NGS <- paste0("TE_df_ls_NGS_",Site_ID,".rds")
  syc_metrics_NGS <- cal_syc_metrics_all_pairs(file_name_NGS)
  # Combine these metrics together
  syc_metrics_site <- as.data.frame(rbind(syc_metrics_full_TS,
                                          syc_metrics_GS,
                                          syc_metrics_NGS))
  # Add Site ID
  syc_metrics_site$Site_ID <- Site_ID
  # Move Site_ID to the first column
  syc_metrics_site <- syc_metrics_site[, c("Site_ID", setdiff(names(syc_metrics_site), "Site_ID"))]
  # Add time period: full-TS, GS, or NGS
  syc_metrics_site$GS <- c("FT","GS","NGS")
  rownames(syc_metrics_site) <- NULL
  # Add this site to all sites
  syc_metrics_all_sites <- rbind(syc_metrics_all_sites,syc_metrics_site)
  print(arrayid)
}

# Output this result df
write.csv(syc_metrics_all_sites,"03_PN_construction/Results/Syc_metrics_all_sites.csv")