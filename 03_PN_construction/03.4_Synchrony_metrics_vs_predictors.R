# Author: Zhaozhen Chen
# Date: 2025.10.22

# This code is to summarize synchrony metrics vs perdictors

# -------- Global ----------
library(dplyr)
library(tidyr)
library(scales)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Predictor df
predictor_df <- read.csv("00_Data/perdictor_df.csv")

# Source functions
source("05_Visualization/Plotting_functions.R")
source("03_PN_construction/Synchrony_metrics_functions_v2.R")

# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")
# All seasons
seasons <- c("FT","GS","NGS")

# All variable combinations
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# Response variable list
res_var_ls <- c("daily_p_TE","daily_agg_TE","best_lag","mem1","mem2","mem3","mem4","mem5")
# Predictor variable list
pre_var_ls <- c("AI","CH","RD","TSand","elevation","porosity")

# Output path
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_vs_predictors/"

# ---------- Main ------------
# Initialize a list to store all results
results_ls <- list()

# Loop over source and sink pairs
for(pair_id in 1:12){
  source_name <- var_comb$from[pair_id]
  sink_name <- var_comb$to[pair_id]
  
  # Read in synchorny metrics df
  Syc_df <- read.csv(paste0(Syc_metrics_path,"Syc_metrics_df_",source_name,"_",sink_name,".csv")) %>%
    dplyr::select(- X)
  
  # Merge with predictors
  Syc_df_merge <- Syc_df %>%
    # Merge with predictors
    left_join(predictor_df %>% dplyr::select(-X),by="site_ID") %>%
    # Merge with site info
    left_join(Site_info %>% dplyr::select(-X),by=c("site_ID" = "site_id"))
  
  # Loop over response and predictors
  for(res_varname in res_var_ls){
    for(pre_varname in pre_var_ls){
      results_ls[[length(results_ls)+1]] <- summarize_syc_stats(Syc_df_merge,source_name,sink_name,res_varname,pre_varname)
    }
  }
  message(pair_id)
}

# Combine all results into a df
results_all <- bind_rows(results_ls)
results_all <- results_all %>%
  mutate(pair = paste0(source,"→",sink)) %>%
  # Rename predictors
  mutate(predictor = if_else(predictor == "elevation","ELE",predictor),
         predictor = if_else(predictor == "porosity","phi",predictor))

# Output this df
write.csv(results_all,paste0(Output_path,"Regression_summary.csv"))

# Explore the results -------------------
# Heat map of adjusted R2, grouped by 12 pairs of source and sink combinations
# Replace negative adjusted R2 with 0
results_all <- results_all %>%
  mutate(R2_adj = if_else(R2_adj > 0,R2_adj,0))

for(season in c("FT","GS","NGS")){
  # Make heatmaps
  g_heatmap_R2 <- heatmap_plot("R2_adj",season,"YlOrRd",direction = 1)
  g_heatmap_k <- heatmap_plot("k",season,"RdBu",direction = -1)
  # Combine the two maps
  g_heatmaps <- plot_grid(g_heatmap_R2,g_heatmap_k,nrow=1,
                          align = "hv")
  print_g(g_heatmaps,paste0("Heatmaps_",season),18,10)
}








  
  


