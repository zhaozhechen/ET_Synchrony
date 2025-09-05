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
Output_path <- "03_PN_construction/Results/"

#season_color <- brewer.pal(3,"Set2")

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
         GS = as.factor(GS)) %>%
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

varname <- "p_TE_psi_to_ET"
palette_name <- "YlOrRd"
var_title <- "Peak daily TE (%)"

# Make maps of the target variable for FT, GS, and NGS
map_3 <- season3_syc_map(Syc_metrics_df,varname,palette_name)










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
# Peak TE
g_p_TE_psi_vs_VPD <- scatter_vars(df = Syc_metrics_df,
                                  "p_TE_psi_to_ET","p_TE_VPD_to_ET",
                                  "GS",
                                  xtitle = bquote(TEmax~"("~Delta~psi~"\u2192"~Delta~ET~")"),
                                  ytitle = bquote(TEmax~"("~Delta~VPD~"\u2192"~Delta~ET~")"),
                                  my_color = season_color)
g_p_TE_psi_vs_TA <- scatter_vars(df = Syc_metrics_df,
                                 "p_TE_psi_to_ET","p_TE_TA_to_ET",
                                 "GS",
                                 xtitle = bquote(TEmax~"("~Delta~psi~"\u2192"~Delta~ET~")"),
                                 ytitle = bquote(TEmax~"("~Delta~'T'[air]~"\u2192"~Delta~ET~")"),
                                 my_color = season_color)
g_p_TE_VPD_vs_TA <- scatter_vars(df = Syc_metrics_df,
                                 "p_TE_VPD_to_ET","p_TE_TA_to_ET",
                                 "GS",
                                 xtitle = bquote(TEmax~"("~Delta~VPD~"\u2192"~Delta~ET~")"),
                                 ytitle = bquote(TEmax~"("~Delta~'T'[air]~"\u2192"~Delta~ET~")"),
                                 my_color = season_color)
# Lag
g_p_lag_psi_vs_VPD <- scatter_vars(df = Syc_metrics_df,
                                   "p_lag_psi_to_ET","p_lag_VPD_to_ET",
                                   "GS",
                                   xtitle = bquote(Lag~"("~Delta~psi~"\u2192"~Delta~ET~")"),
                                   ytitle = bquote(Lag~"("~Delta~VPD~"\u2192"~Delta~ET~")"),
                                   my_color = season_color)
g_p_lag_psi_vs_TA <- scatter_vars(df = Syc_metrics_df,
                                  "p_lag_psi_to_ET","p_lag_TA_to_ET",
                                  "GS",
                                  xtitle = bquote(Lag~"("~Delta~psi~"\u2192"~Delta~ET~")"),
                                  ytitle = bquote(Lag~"("~Delta~'T'[air]~"\u2192"~Delta~ET~")"),
                                  my_color = season_color)
g_p_lag_VPD_vs_TA <- scatter_vars(df = Syc_metrics_df,
                                  "p_lag_VPD_to_ET","p_lag_TA_to_ET",
                                  "GS",
                                  xtitle = bquote(Lag~"("~Delta~VPD~"\u2192"~Delta~ET~")"),
                                  ytitle = bquote(Lag~"("~Delta~'T'[air]~"\u2192"~Delta~ET~")"),
                                  my_color = season_color)
# Memory
g_mem_psi_vs_VPD <- scatter_vars(df = Syc_metrics_df,
                                 "mem_psi_to_ET","mem_VPD_to_ET",
                                 "GS",
                                 xtitle = bquote(Memory~"("~Delta~psi~"\u2192"~Delta~ET~")"),
                                 ytitle = bquote(Memory~"("~Delta~VPD~"\u2192"~Delta~ET~")"),
                                 my_color = season_color)
g_mem_psi_vs_TA <- scatter_vars(df = Syc_metrics_df,
                                "mem_psi_to_ET","mem_TA_to_ET",
                                "GS",
                                xtitle = bquote(Memory~"("~Delta~psi~"\u2192"~Delta~ET~")"),
                                ytitle = bquote(Memory~"("~Delta~'T'[air]~"\u2192"~Delta~ET~")"),
                                my_color = season_color)
g_mem_VPD_vs_TA <- scatter_vars(df = Syc_metrics_df,
                                "mem_VPD_to_ET","mem_TA_to_ET",
                                "GS",
                                xtitle = bquote(Memory~"("~Delta~VPD~"\u2192"~Delta~ET~")"),
                                ytitle = bquote(Memory~"("~Delta~'T'[air]~"\u2192"~Delta~ET~")"),
                                my_color = season_color)
# Combine these plots
g_syc_var <- plot_grid(g_p_TE_psi_vs_VPD,g_p_TE_psi_vs_TA,g_p_TE_VPD_vs_TA,
                       g_p_lag_psi_vs_VPD,g_p_lag_psi_vs_TA,g_p_lag_VPD_vs_TA,
                       g_mem_psi_vs_VPD,g_mem_psi_vs_TA,g_mem_VPD_vs_TA,
                       nrow = 3,
                       align = "hv",
                       labels = "auto")
print_g(g_syc_var,"Syc_metrics_variables",
        10,10)


