# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.8.18

# This code is to explore and analyze synchrony metrics

# -------- Global ----------
library(dplyr)
library(tidyr)

# Input path for Synchrony metrics for all sites
Syc_metrics_df <- read.csv("03_PN_construction/Results/Syc_metrics_all_sites.csv")
# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")
# Source plotting functions
source("05_Visualization/Plotting_functions.R")
# Source synchrony functions
source("03_PN_construction/Synchrony_metrics_functions.R")
# Output path for figures
Output_path <- "03_PN_construction/Results/"

season_color <- brewer.pal(3,"Set2")

# --------- Main ---------
Syc_metrics_df <- Syc_metrics_df  %>%
  select(-X) %>%
  # Only keep GS and Non-GS
  filter(GS != "FT") %>%
  # Join by Site_info
  left_join(Site_info %>% select(Site_ID = site_id,IGBP_veg,Soil_type = Description),
            by = "Site_ID") %>%
  mutate(Soil_type = as.factor(Soil_type),
         IGBP_veg = as.factor(IGBP_veg),
         GS = as.factor(GS))

# Comparison of Synchrony metrics from psi, VPD, and TA to ET =================
# peak TE
g_p_TE_psi_to_ET <- Hist_Syc_p_value(Syc_metrics_df,y_varname = "p_TE_psi_to_ET",group_name = "GS",y_title = "TEmax (%)",
                             x_labels = c("GS","Non-GS"),
                             title = bquote(Delta~psi~"\u2192"~Delta~ET),
                             my_color = season_color,y_lim = c(0,20))
g_p_TE_VPD_to_ET <- Hist_Syc_p_value(Syc_metrics_df,y_varname = "p_TE_VPD_to_ET",group_name = "GS",y_title = "TEmax (%)",
                             x_labels = c("GS","Non-GS"),
                             title = bquote(Delta~VPD~"\u2192"~Delta~ET),
                             my_color = season_color,y_lim = c(0,20))
g_p_TE_TA_to_ET <- Hist_Syc_p_value(Syc_metrics_df,y_varname = "p_TE_TA_to_ET",group_name = "GS",y_title = "TEmax (%)",
                            x_labels = c("GS","Non-GS"),
                            title = bquote(Delta~TA~"\u2192"~Delta~ET),
                            my_color = season_color,y_lim = c(0,20))
# Best lag
g_p_lag_psi_to_ET <- Hist_Syc_p_value(Syc_metrics_df,y_varname = "p_lag_psi_to_ET",group_name = "GS",y_title = "Lag (Hours)",
                              x_labels = c("GS","Non-GS"),
                              title = bquote(Delta~psi~"\u2192"~Delta~ET),
                              my_color = season_color,y_lim = c(0,24))
g_p_lag_VPD_to_ET <- Hist_Syc_p_value(Syc_metrics_df,y_varname = "p_lag_VPD_to_ET",group_name = "GS",y_title = "Lag (Hours)",
                              x_labels = c("GS","Non-GS"),
                              title = bquote(Delta~VPD~"\u2192"~Delta~ET),
                              my_color = season_color,y_lim = c(0,24))
g_p_lag_TA_to_ET <- Hist_Syc_p_value(Syc_metrics_df,y_varname = "p_lag_TA_to_ET",group_name = "GS",y_title = "Lag (Hours)",
                             x_labels = c("GS","Non-GS"),
                             title = bquote(Delta~TA~"\u2192"~Delta~ET),
                             my_color = season_color,y_lim = c(0,24))
# Memory
g_mem_psi_to_ET <- Hist_Syc_p_value(Syc_metrics_df,y_varname = "mem_psi_to_ET",group_name = "GS",y_title = "Memory (Hours)",
                            x_labels = c("GS","Non-GS"),
                            title = bquote(Delta~psi~"\u2192"~Delta~ET),
                            my_color = season_color,y_lim = c(0,72))
g_mem_VPD_to_ET <- Hist_Syc_p_value(Syc_metrics_df,y_varname = "mem_VPD_to_ET",group_name = "GS",y_title = "Memory (Hours)",
                            x_labels = c("GS","Non-GS"),
                            title = bquote(Delta~VPD~"\u2192"~Delta~ET),
                            my_color = season_color,y_lim = c(0,72))
g_mem_TA_to_ET <- Hist_Syc_p_value(Syc_metrics_df,y_varname = "mem_TA_to_ET",group_name = "GS",y_title = "Memory (Hours)",
                           x_labels = c("GS","Non-GS"),
                           title = bquote(Delta~TA~"\u2192"~Delta~ET),
                           my_color = season_color,y_lim = c(0,72))

# Combine all these three plots
g_syc <- plot_grid(g_p_TE_psi_to_ET,g_p_TE_VPD_to_ET,g_p_TE_TA_to_ET,
                    g_p_lag_psi_to_ET,g_p_lag_VPD_to_ET,g_p_lag_TA_to_ET,
                    g_mem_psi_to_ET,g_mem_VPD_to_ET,g_mem_TA_to_ET,
                    nrow = 3,
                    align = "hv",
                    labels = "auto")
print_g(g_syc,"Syc_metrics_season",
        10,10)

# Compare synchrony metrics between variables pairs ============================
# This function is to make scatter plot and compare two variables
varname1 <- "p_TE_psi_to_ET"
varname2 <- "p_TE_VPD_to_ET"
df <- Syc_metrics_df
xtitle <- bquote(TEmax~"("~Delta~psi~"\u2192"~Delta~ET~")")
ytitle <- bquote(TEmax~"("~Delta~VPD~"\u2192"~Delta~ET~")")
group_name <- "GS"










# Plot p_TE_psi_to_ET vs p_TE_VPD_to_ET, color coded by something
