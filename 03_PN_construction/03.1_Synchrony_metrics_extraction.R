# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.9.22

# This code is to extract Synchrony metrics (peak TE, memory, lag, aggTE)
# Peak daily TE: maximum significant TE normalized by Shannon entropy of the sink within the first 24 hours
# Memory: the persistence of influence after the maximum strength of synchrony, defined as the interval from the lag of maximum TE to the first subsequent lag where TE falls below the significance threshold.
# Lag: best lag corresponding to peak TE
# Daily aggregate TE: The total magnitude/contribution of predictive influence distributed across 24 hours, defined as the sum of TE for the set of lags with significant TE within 24 hours.
# Also output normalized TE vs lag + Lomb-Scargle periodogram for all sites, for Full TS, GS, and NGS
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
# Output path for TE+LSP plots
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Lag_plots_allyear_3mem/"
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



# Test below ============
Site_ID <- Site_IDs[1]
type <- "full_TS"

i<-1


# Get file name for the TE_df_ls
file_name <- paste0("TE_df_ls_",type,"_",Site_ID,".rds")
# Read in TE_df_ls
TE_df_ls <- readRDS(paste0(TE_df_path,file_name))



# Source and sink name
source_name <- as.character(var_comb$from[i])
sink_name <- as.character(var_comb$to[i])
# Get TE_df for this pair
TE_df_name <- paste0(source_name,"_to_",sink_name)
TE_df <- TE_df_ls[[TE_df_name]]

cal_syc_metrics(TE_df)



# ====================





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