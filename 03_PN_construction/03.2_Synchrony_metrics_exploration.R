# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.9.4

# This code is to explore and analyze synchrony metrics

# -------- Global ----------
library(dplyr)
library(tidyr)
library(sf)

# Input path for Synchrony metrics for all sites
Syc_metrics_df <- read.csv("03_PN_construction/Results/Syc_metrics_all_sites.csv")
# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")
# Input path for Aridity index
Site_info_AI <- read.csv("00_Data/Site_summary.csv")
# Source plotting functions
source("05_Visualization/Plotting_functions.R")
# Source synchrony functions
source("03_PN_construction/Synchrony_metrics_functions.R")

# Make CONUS boundary
# Whole US map
CONUS <- st_read("00_Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]
# Output path for figures
Output_path <- "03_PN_construction/Results/Syc_metrics_summary/"

# Colors for the three seasons
season_color <- brewer.pal(3,"Set2")
# Palette for making maps
palette_name <- "YlOrRd"

# --------- Main ---------
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












