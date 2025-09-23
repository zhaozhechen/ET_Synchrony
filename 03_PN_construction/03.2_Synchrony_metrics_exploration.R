# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.9.22

# This code is to explore and analyze synchrony metrics

# -------- Global ----------
library(dplyr)
library(tidyr)
library(sf)
library(GGally)
library(patchwork)

# Input path for Synchrony metrics for all sites
Syc_metrics_df <- read.csv("03_PN_construction/Results/Syc_metrics_all_sites_5mem.csv")
# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")
# Input path for Aridity index
#Site_info_AI <- read.csv("00_Data/Site_summary.csv")
# Source plotting functions
source("05_Visualization/Plotting_functions.R")
# Source synchrony metrics functions
source("03_PN_construction/Synchrony_metrics_functions_v2.R")

# Make CONUS boundary
# Whole US map
CONUS <- st_read("00_Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]
# Output path for figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/"

# Colors for the three seasons
season_color <- brewer.pal(3,"Set2")
# Palette for making maps
palette_name <- "YlOrRd"

# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")
# All seasons
seasons <- c("FT","GS","NGS")

# --------- Main ---------
Syc_metrics_df <- Syc_metrics_df  %>%
  select(-X) 

# All variable combinations
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

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















if(FALSE){
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













