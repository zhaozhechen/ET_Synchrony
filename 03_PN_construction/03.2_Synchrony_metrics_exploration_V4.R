# Author: Zhaozhe Chen
# Update Date: 2026.3.31

# This code is to explore lag among pairs

# -------- Global ----------
library(dplyr)
library(tidyr)
library(ggpubr)
library(FSA) # For Dunn test
library(forcats)
library(rstatix)
library(terra)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Predictor df, with updated AI_gridded
predictor_df <- read.csv("00_Data/perdictor_df_updated.csv")

# Aridity index map
AI_raster <- rast("00_Data/2022_CONUS_AI.tif")

# Source-sink pairs to test
source_sink_pairs <- read.csv("00_Data/Source-sink_pairs.csv")

# Source plotting functions
source("05_Visualization/Plotting_functions.R")
source("03_PN_construction/Synchrony_metrics_functions_v2.R")

# Output path of figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Source-sink_pair-comparisons_V5/"

season <- "GS"

# Name of target synchrony metrics
res_varname <- "best_lag"
# Title of these metrics
res_title <- "Lag (h)"

# This is the color pallet for paired comparisons
pair_color <- RColorBrewer::brewer.pal(10,"Paired")[c(1,2,7,8,9,10)]
# Colors for scatter and box plots
syc_colors <- pair_color

# ----- Main -------
# Data preprocessing ==================
# Get full syc_metric name
res_name <- paste0(season,"_",res_varname)

# Read in synchrony metric df for target pairs
syc_df1 <- read_syc_df("psi","ET",Syc_metrics_path)
syc_df2 <- read_syc_df("ET","psi",Syc_metrics_path)
syc_df3 <- read_syc_df("VPD","ET",Syc_metrics_path)
syc_df4 <- read_syc_df("ET","VPD",Syc_metrics_path)
syc_df5 <- read_syc_df("TA","ET",Syc_metrics_path)
syc_df6 <- read_syc_df("ET","TA",Syc_metrics_path)
# Combine these df
syc_df <- rbind(syc_df1,syc_df2,syc_df3,syc_df4,syc_df5,syc_df6) %>%
  merge(predictor_df %>%
          rename(site_ID = site_id),by="site_ID") %>%
  mutate(source_sink = factor(source_sink,levels = c("psi_ET","ET_psi","VPD_ET","ET_VPD","TA_ET","ET_TA"))) %>%
  # Regime, based on Dryness Index = PET/P = 1/AI
  # PET/P < 1 (AI>1): Energy-limited
  # PET/P > 1 (AI<1): Water-limited
  mutate(Regime = case_when(AI_gridded < 1 ~ "Water-limited",
                            AI_gridded >= 1 ~ "Energy-limited")) %>%
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
                                       TRUE ~ Koppen_clim_class))%>%
  # Add empty levels to Koppen climate classes for plotting
  mutate(
    Koppen_clim_class = factor(Koppen_clim_class,levels = c("BSh","BSk","BWh","BWk","",
                                                            "Csa","Csb"," ",
                                                            "Cfa","Cwa","  ",
                                                            "Dfa","Dfb","Dsa","Dsb","   ",
                                                            "Dfc"))) %>%
  # Aggregate some IGBP classes, due to small sample size in some of the original classes
  mutate(IGBP_aggregate = case_when(
    IGBP_veg %in% c("ENF","DBF","MF") ~ "Forest",
    IGBP_veg %in% c("CSH","OSH") ~ "Shrubland",
    IGBP_veg %in% c("WSA","SAV","GRA")~"Savanna/Grass",
    IGBP_veg %in% c("CRO","CVM")~ "Cropland"
  )) %>%
  mutate(IGBP_aggregate = factor(IGBP_aggregate,levels = c("Cropland","Savanna/Grass","Shrubland","Forest"))) %>%
  # Add empty levels to IGBP for plotting
  mutate(IGBP_veg = factor(IGBP_veg,levels = c("ENF","DBF","MF","",
                                               "CSH","OSH","  ",
                                               "WSA","SAV","GRA","   ",
                                               "CRO","CVM"))) %>%
  left_join(Site_info %>%
              select(site_ID = site_id,Soil_Type = Description),by="site_ID") %>%
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
  Soil_Group = factor(Soil_Group,levels=c("Fine","Moderately fine","Medium","Moderately coarse","Coarse")))

# Make plot including both lag and peak TE =========================
g_lag_TE_AI <- plot_TE_vs_lag("AI_Class")
g_lag_TE_climate <- plot_TE_vs_lag("Koppen_aggregate")+
  guides(size = "none")
g_lag_TE_IGBP <- plot_TE_vs_lag("IGBP_aggregate")+
  guides(size = "none")
g_lag_TE_soil <- plot_TE_vs_lag("Soil_Group")+
  guides(size = "none")

# Combine these plots
g_lag_TE <- plot_grid(g_lag_TE_AI,g_lag_TE_climate,g_lag_TE_IGBP,g_lag_TE_soil,
                      ncol=1,align="hv")

print_g(g_lag_TE,paste0("Lag_TE_plot_v2"),8,16)

# Distributions of lag time across 6 pairs ===============
g_box <- ggplot(data = syc_df,aes(x=source_sink,y=.data[[res_name]],fill=source_sink))+
  geom_boxplot(outlier.color = "grey")+
  my_theme+
  labs(y = res_title,x="")+
  scale_fill_manual(values = pair_color)+
  scale_x_discrete(labels = c(
    psi_ET = expression(psi %->% ET),
    ET_psi = expression(ET %->% psi),
    VPD_ET = expression(VPD %->% ET),
    ET_VPD = expression(ET %->% VPD),
    TA_ET  = expression(T[air] %->% ET),
    ET_TA  = expression(ET %->% T[air])
  ))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
print_g(g_box,paste0("Box_",res_name),3,4)

# Maps of lag across pairs =====================
# A list of source-sink pairs to plot
sc_pairs <- data.frame(
  source = c("psi","VPD","TA","ET","ET","ET"),
  sink   = c("ET","ET","ET","psi","VPD","TA")
)

# named colors for the 6 pairs
pair_names <- paste0(sc_pairs$source, "_", sc_pairs$sink)
pair_cols  <- syc_colors[seq_len(nrow(sc_pairs))]
names(pair_cols) <- pair_names

# initialize lists
g_map_ls <- list()
syc_df_tmp <- list()
g_scatter_box_ls <- list()
g_scatter_box_SI_ls <- list()
g_scatter_ls <- list()

for (i in 1:nrow(sc_pairs)) {
  source_name <- sc_pairs$source[i]
  sink_name   <- sc_pairs$sink[i]
  pair_name   <- paste0(source_name, "_", sink_name)
  
  # Remove NA
  syc_df_sub <- syc_df %>%
    filter(source_sink == pair_name) %>%
    filter(!is.na(.data[[res_name]]))
  
  syc_df_tmp[[i]] <- syc_df_sub
  
  # Scatter
  g_scatter <- scatter_vars(
    syc_df_sub, "AI_gridded", res_name, "source_sink",
    "Aridity Index", res_title, syc_colors[i]
  ) +
    theme(legend.position = "none")
  
  g_scatter_ls[[i]] <- g_scatter
  
  # Map
  map_color <- RColorBrewer::brewer.pal(n = 11, name = "PuBu")
  
  # pair-specific limits
 # if (pair_name %in% c("psi_ET", "ET_psi")) {
#    color_limits_i <- c(0, 24)
#    n_breaks_i <- 5
#  } else {
#    color_limits_i <- c(0, 0.5)
#    n_breaks_i <- 6
#  }
  
  g_map <- syc_map(
    df = syc_df_sub,
    varname = res_name,
    colors = map_color,
    legend_title = res_title,
    g_title = "",
    color_limits = c(0,24),
    n_breaks = 5,
    end_marks = "right",
    base_raster = AI_raster
  ) +
    theme(legend.position = "bottom")
  
  g_map_ls[[i]] <- g_map
}

# Combine 6 maps ====
g_map <- plot_grid(
  plotlist = g_map_ls,
  nrow = 2,
  align = "hv",
  labels = "auto"
)

print_g(g_map, paste0("Lag_maps_", res_name), 21, 7)



# bind all pairs together for faceted box/violin plots
syc_df_tmp <- dplyr::bind_rows(syc_df_tmp)

# Main text plots: aggregated groups -----------------------
g_box_AI_facet <- plot_box_groups_facet(
  syc_df_tmp, res_name, "AI_Class",
  fill_colors = pair_cols,
  x_title = res_title, y_title = ""
)

g_box_climate_agg_facet <- plot_box_groups_facet(
  syc_df_tmp, res_name, "Koppen_aggregate",
  fill_colors = pair_cols,
  x_title = res_title, y_title = ""
)

g_box_IGBP_agg_facet <- plot_box_groups_facet(
  syc_df_tmp, res_name, "IGBP_aggregate",
  fill_colors = pair_cols,
  x_title = res_title, y_title = ""
)

g_box_soil_agg_facet <- plot_box_groups_facet(
  syc_df_tmp, res_name, "Soil_Group",
  fill_colors = pair_cols,
  x_title = res_title, y_title = ""
)

g_scatter_box_ls <- list(
  g_box_AI_facet,
  g_box_climate_agg_facet,
  g_box_IGBP_agg_facet,
  g_box_soil_agg_facet
)

# SI plots: scatter + detailed groups ----------------------
g_box_climate_facet <- plot_box_groups_facet(
  syc_df_tmp, res_name, "Koppen_clim_class",
  fill_colors = pair_cols,
  x_title = res_title, y_title = "",
  box_violin = "Box"
)

g_box_IGBP_facet <- plot_box_groups_facet(
  syc_df_tmp, res_name, "IGBP_veg",
  fill_colors = pair_cols,
  x_title = res_title, y_title = "",
  box_violin = "Box"
)

g_box_soil_facet <- plot_box_groups_facet(
  syc_df_tmp, res_name, "Soil_Type",
  fill_colors = pair_cols,
  x_title = res_title, y_title = "",
  box_violin = "Box"
)

g_scatter_box_SI_ls <- c(
  list(
    g_box_climate_facet,
    g_box_IGBP_facet,
    g_box_soil_facet
  )
)

# Combine all scatter and box plots (main text)
g_scatter_box <- plot_grid(plotlist = g_scatter_box_ls,
                           ncol=1,align="hv",labels="auto")
print_g(g_scatter_box,"Lag_compare_agg",12,24)

# Combine all scatter and box plots (SI figure)
g_scatter_box_SI <- plot_grid(plotlist = g_scatter_box_SI_ls,
                              ncol=1,align="hv",labels = "auto")
print_g(g_scatter_box_SI,"Lag_compare_SI",12,24)
