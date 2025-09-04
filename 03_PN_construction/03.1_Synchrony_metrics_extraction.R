# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.9.4

# This code is to extract Synchrony metrics (peak TE, memory, lag)
# Peak TE: maximum significant TE normalized by Shannon entropy of the sink
# Memory: time required for TE to become insignificant
# Lag: best lag corresponding to peak TE
# Note: Focus on full TE, GS, and NGS for now. No years results included

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
source("03_PN_construction/Synchrony_metrics_functions.R")
source("05_Visualization/Plotting_functions.R")
# Output path for TE+LSP plots
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/TE_LS_plots/"
my_color <- brewer.pal(5,"Set1")

# Set parallel session
plan(multisession,workers = availableCores()-1)

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

# Process for all sites
with_progress({
  # Initiate a progressor
  p <- progressor(along = Site_IDs)
  syc_metrics_all_sites_ls <- future_lapply(Site_IDs,function(Site_ID){
    result <- syc_site(Site_ID)
    p() # Update progress
    return(result)
  })
})

# Combine all output
syc_metrics_all_sites <- do.call(rbind,syc_metrics_all_sites_ls)

# Output this combined df
write.csv(syc_metrics_all_sites,"03_PN_construction/Results/Syc_metrics_all_sites.csv")