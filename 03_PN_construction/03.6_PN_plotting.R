# Author: Zhaozhe Chen
# Update Date: 2026.4.7

# This code is to plot Process Network Chord diagram
# Reference: https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html#directional-relations

# ----- Global -----------
library(dplyr)
library(circlize) # For chord diagram
library(RColorBrewer)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"
# Updated site info
site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")
# Read predictor df (for AI)
predictor_df <- read.csv("00_Data/perdictor_df_updated.csv")

# Source functions for data processing
source("03_PN_construction/Synchrony_metrics_functions_v2.R")
# Source functions for plotting
source("05_Visualization/Plotting_functions.R")

# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")
# All variable combinations
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# Color for the ribbons 
my_colors <- RColorBrewer::brewer.pal(10,"Paired")[c(3,1,7,9)]
cols <- c(ET=my_colors[1], psi=my_colors[2], VPD=my_colors[3], TA=my_colors[4]) 

# Output path for figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Aggregated_Chord_diagrams_V2/"

# target synchrony metric name to be plotted in the chord diagram
#target_syc_metric_name <- "daily_p_TE"
target_syc_metric_name <- "daily_agg_TE"
# Season to be plotted
season <- "GS"

# ------ Main -------
# Data processing =============================
Site_info <- predictor_df %>%
  # Aridity classes, based on Huang et al. 2016
  mutate(AI_Class = case_when(AI_gridded <0.2 ~ "Arid",
                              (AI_gridded >= 0.2 & AI_gridded < 0.5) ~ "Semiarid",
                              (AI_gridded >= 0.5 & AI_gridded < 0.65) ~ "Semihumid",
                              AI_gridded >= 0.65 ~ "Humid")) %>%
  mutate(AI_Class = factor(AI_Class,levels=c("Arid","Semiarid","Semihumid","Humid")))%>%
  # Aggregate some Koppen climate classes, due to small sample size in some of the original classes
  mutate(
    Koppen_aggregate = case_when(
      Koppen_clim_class %in% c("Bsh","Bsk","Bwh","Bwk") ~ "Dry",
      Koppen_clim_class %in% c("Cfa","Cwa")              ~ "Humid_subtropical",
      Koppen_clim_class %in% c("Csa","Csb")              ~ "Mediterranean",
      Koppen_clim_class %in% c("Dfa","Dfb","Dsa","Dsb")  ~ "Humid_continental",
      Koppen_clim_class == "Dfc"                         ~ "Subarctic",
      TRUE ~ Koppen_clim_class
    )
  ) %>%
  mutate(Koppen_aggregate = factor(Koppen_aggregate,levels = c("Dry","Mediterranean","Humid_subtropical",
                                                               "Humid_continental","Subarctic"))) %>%
  # Update climate classes names
  mutate(Koppen_clim_class = case_when(Koppen_clim_class == "Bsh"~"BSh",
                                       Koppen_clim_class == "Bsk"~"BSk",
                                       Koppen_clim_class == "Bwh"~"BWh",
                                       Koppen_clim_class == "Bwk"~"BWk",
                                       TRUE ~ Koppen_clim_class)) %>%
  # Aggregate some IGBP classes, due to small sample size in some of the original classes
  mutate(IGBP_aggregate = case_when(
    IGBP_veg %in% c("ENF","DBF","MF") ~ "Forest",
    IGBP_veg %in% c("CSH","OSH") ~ "Shrubland",
    IGBP_veg %in% c("WSA","SAV","GRA")~"Savanna/Grass",
    IGBP_veg %in% c("CRO","CVM")~ "Cropland"
  )) %>%
  mutate(IGBP_aggregate = factor(IGBP_aggregate,levels = c("Cropland","Savanna/Grass","Shrubland","Forest"))) %>%
  left_join(site_info %>%
              select(site_id,Soil_Type = Description),by="site_id") %>%
  mutate(Soil_Type = factor(Soil_Type,levels=rev(c("Sand","Loamy sand",
                                                   "    ","Sandy loam",
                                                   "   ","Loam","Silt loam",
                                                   "  ","Clay loam","Silty clay loam",
                                                   " ","Silty clay","Clay")))) %>%
  # Ref: USDA The Soil Survey Manual 2017 p123-126
  mutate(Soil_Group = case_when(
    Soil_Type %in% c("Clay","Silty clay") ~ "Fine",
    Soil_Type %in% c("Silty clay loam","Clay loam") ~ "Moderately fine",
    Soil_Type %in% c("Loam","Silt loam") ~ "Medium",
    Soil_Type %in% c("Sandy loam") ~ "Moderately coarse",
    Soil_Type %in% c("Loamy sand","Sand") ~ "Coarse",
    TRUE ~ NA_character_
  ),
  Soil_Group = factor(Soil_Group,levels=c("Fine","Moderately fine","Medium","Moderately coarse","Coarse"))) %>%
  mutate(CONUS = "CONUS")

# Across CONUS ==============
make_group_chord_diagram(group_col = "CONUS",
                         title = "CONUS",
                         panel_w = 4,
                         panel_h = 4,
                         Site_info = Site_info,
                         var_comb = var_comb,
                         Syc_metrics_path = Syc_metrics_path,
                         season = season,
                         target_syc_metric_name = target_syc_metric_name,
                         Output_path = Output_path,
                         cols = cols,
                         ncol_fixed = 1,
                         var_order = c("ET","psi","VPD","TA"),
                         sector_label_cex = 1.5)

# Group the sites by IGBP vegetation type =======================
make_group_chord_diagram(group_col = "IGBP_aggregate",
                         title = "IGBP",
                         Site_info = Site_info,
                         var_comb = var_comb,
                         Syc_metrics_path = Syc_metrics_path,
                         season = season,
                         target_syc_metric_name = target_syc_metric_name,
                         Output_path = Output_path,
                         cols = cols,
                         ncol_fixed = 4,
                         var_order = c("ET","psi","VPD","TA"))

# Group the sites by Koppen climate classification =======================
make_group_chord_diagram(group_col = "Koppen_aggregate",
                         title = "Koppen_climate",
                         Site_info = Site_info,
                         var_comb = var_comb,
                         Syc_metrics_path = Syc_metrics_path,
                         season = season,
                         target_syc_metric_name = target_syc_metric_name,
                         Output_path = Output_path,
                         cols = cols,
                         var_order = c("ET","psi","VPD","TA"))

# Group the sites by Soil type =======================
make_group_chord_diagram(group_col = "Soil_Group",
                         title = "Soil_texture",
                         Site_info = Site_info,
                         var_comb = var_comb,
                         Syc_metrics_path = Syc_metrics_path,
                         season = season,
                         target_syc_metric_name = target_syc_metric_name,
                         Output_path = Output_path,
                         cols = cols,
                         var_order = c("ET","psi","VPD","TA"))

# Group the sites by AI class =======================
make_group_chord_diagram(group_col = "AI_Class",
                         title = "AI_class",
                         Site_info = Site_info,
                         var_comb = var_comb,
                         Syc_metrics_path = Syc_metrics_path,
                         season = season,
                         target_syc_metric_name = target_syc_metric_name,
                         Output_path = Output_path,
                         cols = cols,
                         ncol_fixed = 4,
                         var_order = c("ET","psi","VPD","TA"))

# For each site =====================
make_group_chord_diagram(group_col = "site_id",
                         panel_w = 4,
                         panel_h = 4,
                         title = "Sites",
                         Site_info = Site_info,
                         var_comb = var_comb,
                         Syc_metrics_path = Syc_metrics_path,
                         season = season,
                         target_syc_metric_name = target_syc_metric_name,
                         Output_path = Output_path,
                         cols = cols,
                         var_order = c("ET","psi","VPD","TA"),
                         sector_label_cex = 2)


