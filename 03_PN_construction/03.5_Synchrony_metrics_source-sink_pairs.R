# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.11.4

# This code is to plot heatmap of delta_response between source-sink pairs vs predictors

# ----- Gloabl --------
library(dplyr)
library(tidyr)
library(scales)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Predictor df
predictor_df <- read.csv("00_Data/perdictor_df_updated.csv")

# Source-sink pairs to test
source_sink_pairs <- read.csv("00_Data/Source-sink_pairs.csv")

# Source plotting functions
source("05_Visualization/Plotting_functions.R")
source("03_PN_construction/Synchrony_metrics_functions_v2.R")

# Response variable list
res_var_ls <- c("daily_p_TE","daily_agg_TE","best_lag","mem1","mem2","mem3","mem4","mem5")
# Predictor variable list
pre_var_ls <- c("AI","CH_GLAD","RD","AVG_CLAY","AVG_SILT","elevation","porosity")

Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Source-sink_pair-comparisons_updated/"

# ---- Main -------
# Initialize a list to store all regression results
results_delta_ls <- list()

for(res_varname in res_var_ls){
  for(arrayid in 1:nrow(source_sink_pairs)){
    source_name1 <- source_sink_pairs$Source1[arrayid]
    sink_name1 <- source_sink_pairs$Sink1[arrayid]
    source_name2 <- source_sink_pairs$Source2[arrayid]
    sink_name2 <- source_sink_pairs$Sink2[arrayid]
    
    # Make merged syc metrics, only keep target response variable (e.g., daily_p_TE), and calculate delta_response
    pair_syc_df_merged <- merge_pair_syc_df(source_name1,sink_name1,source_name2,sink_name2,res_varname)
    
    # If delta in lags between two source-sink pairs is < 0, then it should be plus 24
    if(res_varname == "best_lag"){
      pair_syc_df_merged <- pair_syc_df_merged %>%
        mutate(delta_response = if_else(delta_response < 0,delta_response + 24,delta_response))
    }
    
    # Get regression results
    results_delta_ls[[length(results_delta_ls) + 1]] <- 
      summarize_delta_stats(pair_syc_df_merged,source_name1,sink_name1,source_name2,sink_name2,res_varname,pre_var_ls)
  }
}
# Combine all results
results_delta_all <- bind_rows(results_delta_ls) %>%
  # Rename predictors
  mutate(predictor = if_else(predictor == "elevation","ELE",predictor),
         predictor = if_else(predictor == "porosity","phi",predictor)) %>%
  mutate(R2_adj = if_else(R2_adj > 0,R2_adj,0))

# Make a heat maps
for(season in c("FT","GS","NGS")){
  g_heatmap_R2 <- plot_delta_heatmaps(results_delta_all,"R2_adj",season,"YlOrRd",direction = 1)
  g_heatmap_k <- plot_delta_heatmaps(results_delta_all,"k",season,"RdBu",direction = -1)
  # Combine the two maps
  g_heatmaps <- plot_grid(g_heatmap_R2,g_heatmap_k,nrow=2,
                          align = "hv")
  print_g(g_heatmaps,paste0("Heatmaps_",season),16,14)
}




