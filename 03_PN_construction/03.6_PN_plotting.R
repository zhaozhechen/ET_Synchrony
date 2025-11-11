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

# node order around the circle
var_order <- c("ET","psi","VPD","TA")
# Color for the ribbons 
my_colors <- RColorBrewer::brewer.pal(4,"Set3")
cols <- c(ET=my_colors[1], psi=my_colors[2], VPD=my_colors[4], TA=my_colors[3]) 

# Only focus on growing season
#season <- "GS"
# This is the variable name for the target synchrony metric to be plotted in the chord diagram
#target_syc_metric_name <- "daily_p_TE"

# Output path for figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Aggregated_Chord_diagrams/"

# ------ Main -------

# This is the index for Site ID to be included in this aggregation
Site_ID_index <- c(1,2,3)
Site_ls <- Site_info$site_id[Site_ID_index]

# Get df of mean synchrony metrics across target sites
plot_df <- get_target_syc_metric(Site_ls,var_comb,Syc_metrics_path,season = "GS",target_syc_metric_name = "daily_p_TE")

# Make chord diagram
plot_chord_diagram(plot_df,cols = cols,var_order=var_order)


figure_filename <- "Test"

w=2200
h=2000
r=200
pdf(paste0(Output_path,"/",figure_filename,".pdf"),
    width=w,height=h)





