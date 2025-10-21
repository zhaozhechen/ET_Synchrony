# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.10.21

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

# Input path for TE_df
TE_df_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/TE_df/"
# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")
# Source functions
source("03_PN_construction/Synchrony_metrics_functions_v2.R")
source("05_Visualization/Plotting_functions.R")

# Output path for syc metrics
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

max_lag <- 72 # Maximum lag to consider
# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")

# All variable combinations (12 in total)
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# Sites to process
Site_IDs <- Site_info$site_id
# Types to process
type_ls <- c("full_TS","GS","NGS")
# The corresponding type names
typename_ls <- c("FT","GS","NGS")

# ------- Main ------
# Determine which pair of variable to process
for(pair_id in 1:12){
  # Source and sink name
  source_name <- as.character(var_comb$from[pair_id])
  sink_name <- as.character(var_comb$to[pair_id])
  
  # Initialize a df to store all output syc metrics for one pair of variables, across all sites
  syc_metrics_all <- data.frame()
  # Loop over all sites
  for(arrayid in 1:length(Site_IDs)){
    Site_ID <- Site_IDs[arrayid]
    
    # Initialize a vector to store syc metrics at this site, across three types
    syc_metrics_3types <- c()
    # Loop over three types
    for(type_id in 1:3){
      # Determine which type (full TS, GS, or NGS) to process
      type <- type_ls[type_id]
      # The corresponding type name
      typename <- typename_ls[type_id]
      
      # Get file name for the TE_df_ls
      file_name <- paste0("TE_df_ls_",type,"_",Site_ID,".rds")
      # Read in TE_df_ls
      TE_df_ls <- readRDS(paste0(TE_df_path,file_name))
      # Get TE_df for this pair
      TE_df_name <- paste0(source_name,"_to_",sink_name)
      TE_df <- TE_df_ls[[TE_df_name]]
      
      # Check if this TE_df is valid
      if(nrow(TE_df) < max_lag){
        # All NA
        syc_metrics <- rep(NA,16)
        names(syc_metrics) <- c(
          "daily_p_TE","mean_TE","mean_sig_TE","daily_agg_TE","var_TE",
          "sk_TE","best_lag","mem1","mem2","mem3","mem4","mem5",
          "cv_TE","ac1_TE","RB_id_TE","H_TE")
      }else{
        # Use constant Hy values for all lags so that TEcrit after converting to % doesn't vary across lags
        TE_df$Hy <- TE_df$Hy[1]
        # Calculate Synchrony metrics if valid
        syc_metrics <- cal_syc_metrics(TE_df)
      }
      
      # Prefix the syc metric names with typename
      names(syc_metrics) <- paste0(typename,"_",names(syc_metrics))
      
      # Combine 3 types into 1 row
      syc_metrics_3types <- c(syc_metrics_3types,syc_metrics)
    }
    
    # Combine syc metrics for all sites
    syc_metrics_all <- rbind(syc_metrics_all,
                             data.frame(site_ID = Site_ID,
                                        t(syc_metrics_3types)))
  }
  
  # Output all syc metrics for this pair of variables
  write.csv(syc_metrics_all,paste0(Output_path,"Syc_metrics_df_",source_name,"_",sink_name,".csv"))
  message("Complete Pair",pair_id)
}



