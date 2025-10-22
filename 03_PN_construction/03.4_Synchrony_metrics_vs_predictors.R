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


# ---------- Main ------------

pair_id <- 1
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

res_varname <- res_var_ls[1]
pre_varname <- pre_var_ls[1]






  
  


