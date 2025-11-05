# Author: Zhaozhe Chen
# Date: 2025.11.4

# This code is to plot Process Network Chord diagram

# ----- Global -----------
library(dplyr)
library(circlize) # For chord diagram

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")
# All seasons
seasons <- c("FT","GS","NGS")

# All variable combinations
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# ------ Main -------
Site_ID <- "US-A32"
season <- "GS"
# This is the variable name for the target synchrony metric
target_syc_metric_name <- "daily_p_TE"

# Get target synchrony metric between all source-sink pairs, for the target season
for(arrayid in 1:nrow(var_comb)){
  # Source variable name
  source_name <- var_comb$from[arrayid]
  # Sink variable name
  sink_name <- var_comb$to[arrayid]
  
  # Read in synchrony metrics df
  syc_df <- read.csv(paste0(Syc_metrics_path,"Syc_metrics_df_",source_name,"_",sink_name,".csv"))
  
  # Extract only the target variable
  
}

