# Author: Zhaozhe Chen
# Date: 2025.11.11

# This code is to plot Process Network Chord diagram
# Reference: https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html#directional-relations

# ----- Global -----------
library(dplyr)
library(circlize) # For chord diagram
library(RColorBrewer)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"
# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

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
my_colors <- RColorBrewer::brewer.pal(4,"Set3")
cols <- c(ET=my_colors[1], psi=my_colors[2], VPD=my_colors[4], TA=my_colors[3]) 

# Output path for figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Aggregated_Chord_diagrams/"

# ------ Main -------
# Add Climate zone to sites that have NA Koppen_clim_class
Site_info$Koppen_clim_class[Site_info$site_id == "US-Akn"] <- "Cfa"
Site_info$Koppen_clim_class[Site_info$site_id == "US-Cop"] <- "Bsk"
Site_info$Koppen_clim_class[Site_info$site_id == "US-Lin"] <- "Csa"

# Group the sites by IGBP vegetation type =======================
make_group_chord_diagram(group_col = "IGBP_veg",
                         title = "IGBP",
                         Site_info = Site_info,
                         var_comb = var_comb,
                         Syc_metrics_path = Syc_metrics_path,
                         season = "GS",
                         target_syc_metric_name = "daily_p_TE",
                         Output_path = Output_path,
                         cols = cols,
                         var_order = c("ET","psi","VPD","TA"))

# Group the sites by Koppen climate classification =======================
make_group_chord_diagram(group_col = "Koppen_clim_class",
                         title = "Koppen_climate",
                         Site_info = Site_info,
                         var_comb = var_comb,
                         Syc_metrics_path = Syc_metrics_path,
                         season = "GS",
                         target_syc_metric_name = "daily_p_TE",
                         Output_path = Output_path,
                         cols = cols,
                         var_order = c("ET","psi","VPD","TA"))





















