# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.10.22

# This code is to explore and analyze synchrony metrics

# -------- Global ----------
library(dplyr)
library(tidyr)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Predictor df
predictor_df <- read.csv("00_Data/perdictor_df.csv")

# Source plotting functions
source("05_Visualization/Plotting_functions.R")

# Output path for figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_vs_predictors/"

# Colors for the three seasons
season_color <- brewer.pal(3,"Set2")
# Palette for making maps
palette_name <- "YlOrRd"

# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")
# All seasons
seasons <- c("FT","GS","NGS")

# All variable combinations
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# Determine whether to make plots
If_plot <- FALSE

# --------- Main ---------
# Response variable list
res_var_ls <- c("daily_p_TE","daily_agg_TE","best_lag","mem1","mem2","mem3","mem4","mem5")
# Predictor variable list
pre_var_ls <- c("AI","CH","RD","TSand","elevation","porosity")

for(pre_id in 1:length(pre_var_ls)){
  # Variable name for predictor
  pre_varname <- pre_var_ls[pre_id]
  # Loop over each pair of source and sink variables
  # Initialize a list to store the figures
  g_ls <- list()
  for(pair_id in 1:12){
    # Determine which variable pair to process
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
    
    # Loop over the response variables
    for(res_id in 1:length(res_var_ls)){
      res_varname <- res_var_ls[res_id]
    
      # Make a figure
      g_scatter <- syc_scatter(Syc_df_merge,res_varname,pre_varname)+
        ggtitle(paste(source_name,"→",sink_name))
      # Store this figure in the figure list
      g_ls[[length(g_ls)+1]] <- g_scatter
    }
  }
  
  # Combine all figures for this predictor into 1 plot
  g_all <- plot_grid(plotlist = g_ls,ncol = length(res_var_ls),
                     align = "hv")
  # Output this figure
  print_g(g_all,paste0("Syc_metrics_vs_",pre_varname),
          24,36)
  message(pre_id,"out of",length(pre_var_ls))
}


# Previous plotting codes =======================================================
if(FALSE){
  # Compare 5 memory options ========
  # Loop each season and variable pair
  for(i in 1:nrow(var_comb)){
    source_name <- var_comb$from[i]
    sink_name <- var_comb$to[i]
    # Get memory columns names
    mem_names <- paste0(c("mem1","mem2","mem3","mem4","mem5"),"_",source_name,"_to_",sink_name)
    # Initiliate a list to store three seasons
    g_seasons <- list()
    for(season in seasons){
      # Get a subset of these memory columns for this season
      df_tmp <- Syc_metrics_df %>%
        filter(GS == season) %>%
        select(Site_ID,all_of(mem_names))
      # Rename the columns
      names(df_tmp) <- c("Site_ID","Mem1","Mem2","Mem3","Mem4","Mem5")
      # Pivot to long data
      df_tmp_long <- df_tmp %>%
        pivot_longer(cols = starts_with("Mem"),
                     names_to = "Metric",
                     values_to = "value")
      # Create scatter plot matrix
      g_matrix <- ggpairs(df_tmp[,c(2:6)])+
        ggtitle(paste0("Memory comparison (",source_name," → ",sink_name,"; ",season,")"))+
        my_theme
      g_seasons[[season]] <- grob_to_ggplot(g_matrix)
    }
    
    
    # Output plots for this variable pair
    g <- wrap_plots(g_seasons,nrow=1)
    print_g(g,paste0("Mem5_comparison/Mem4_",source_name,"_to_",sink_name),
            15,5)
    
    message(i)
  }
  
  # Preprocessing of synchrony df -------------------
  Syc_metrics_df <- Syc_metrics_df  %>%
    select(-X) %>%
    # Join by Site_info
    left_join(Site_info %>% rename(Site_ID = site_id,Soil_type = Description),
              by = "Site_ID") %>%
    # Join by AI
    left_join(Site_info_AI %>% select(Site_ID,AI),by="Site_ID") %>%
    mutate(Soil_type = as.factor(Soil_type),
           IGBP_veg = as.factor(IGBP_veg),
           GS = as.factor(GS),
           Koppen_clim_class = as.factor(Koppen_clim_class)) %>%
    # Classify AI into five gradients
    mutate(
      AI_level = case_when(
        AI < 0.05 ~ "Hyperarid",
        AI >= 0.05 & AI < 0.2 ~ "Arid",
        AI >= 0.2 & AI < 0.5 ~ "Semiarid",
        AI >= 0.5 & AI < 0.65 ~ "Semihumid",
        AI >= 0.65 ~ "Humid"
      ),
      AI_level = factor(AI_level,
                        levels=c("Hyperarid","Arid","Semiarid","Semihumid","Humid"),
                        ordered = TRUE)
    )
  
  # Make plots for target variable ---------------------
  
  # List of syc metrics to plot
  varname_ls <- names(Syc_metrics_df)[2:61]
  
  # Loop over each variable
  for(varname in varname_ls){
    plot_syc_all(Syc_metrics_df,varname,palette_name)  
  }  
}













