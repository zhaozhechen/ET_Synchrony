# Author: Zhaozhe
# Date: 2026.2.26

# This code is to make a map of dominant driver to ET

# -------- Global ----------
library(dplyr)
library(tidyr)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Source plotting functions
source("05_Visualization/Plotting_functions.R")
source("03_PN_construction/Synchrony_metrics_functions_v2.R")

# Output path of figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Dominant_drivers/"

season <- "GS"

# synchrony variable name
syc_varname <- paste0(season,"_daily_p_TE")

# Read in synchrony metric df for target pairs
syc_df1 <- read_syc_df("psi","ET",Syc_metrics_path) %>%
  select(site_ID,all_of(syc_varname),"source_sink")
syc_df2 <- read_syc_df("VPD","ET",Syc_metrics_path) %>%
  select(site_ID,all_of(syc_varname),"source_sink")
syc_df3 <- read_syc_df("TA","ET",Syc_metrics_path) %>%
  select(site_ID,all_of(syc_varname),"source_sink")

# Combine
syc_all <- bind_rows(syc_df1,syc_df2,syc_df3) %>%
  # Only keep driver name (to ET)
  mutate(DriverName = sub("_ET","",source_sink)) %>%
  select(site_ID,DriverName,TE = all_of(syc_varname)) %>%
  pivot_wider(names_from = DriverName,values_from = TE) %>%
  rowwise() %>%
  mutate(
    # vector of values
    vals = list(c(psi = psi, VPD = VPD, TA = TA)),
    
    Driver = if (all(is.na(unlist(vals)))) {
      NA_character_
    } else {
      names(which.max(replace(unlist(vals), is.na(unlist(vals)), -Inf)))
    },
    
    Process = if (all(is.na(unlist(vals)))) {
      NA_character_
    } else {
      paste(
        names(sort(unlist(vals), decreasing = TRUE, na.last = TRUE)),
        collapse = "-"
      )
    }
  ) %>%
  select(-vals) %>%
  ungroup() %>%
  left_join(Site_info %>%
              select(site_ID = site_id,
                     latitude,longitude),
            by="site_ID")
  
my_color <- RColorBrewer::brewer.pal(9,"Set1")
g_Driver <- syc_map_disc(syc_all,"Driver","Driver","")

my_color <- RColorBrewer::brewer.pal(9,"Set1")
g_Process <- syc_map_disc(syc_all,"Process","Process","")

# Output these two figures
g <- plot_grid(g_Driver,g_Process,nrow=2)
print_g(g,"Driver+Process",8,8)
