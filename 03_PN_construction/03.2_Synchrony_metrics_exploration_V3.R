# Author: Zhaozhe Chen
# Date: 2026.5.1

# This code is to analyze Delta (TE(driver->ET) - TE(ET->driver))

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
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS_LAI_filtered.csv")

# Predictor df
predictor_df <- read.csv("00_Data/perdictor_df_updated.csv")

# Aridity index map
AI_raster <- rast("00_Data/2022_CONUS_AI.tif")

# Source plotting functions
source("05_Visualization/Plotting_functions.R")
source("03_PN_construction/Synchrony_metrics_functions_v2.R")

Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Source-sink_pair-comparisons_V4/"

# Season
season <- "GS"

# Determine which syc metric to process
arrayid <- 1
# Name of target synchrony metrics
res_varname <- c("daily_p_TE","daily_agg_TE","best_lag")[arrayid]
# Title of these metrics
res_title <- c("Peak TE (%)","Daily aggregated TE (%)","Lag (h)")[arrayid]

# This is the color pallet for paired comparisons
pair_color <- RColorBrewer::brewer.pal(10,"Paired")[c(1,2,7,8,9,10)]

# Title for Delta and TE total
Delta_title <- "Directional asymmetry (ΔTE, %)"
TEtotal_title <- "Bidirectional coupling strength (TEtotal, %)"

# ------- Main --------
# Data processing ===============================
# Get full syc_metric name
res_name <- paste0(season,"_",res_varname)

# Source and sink pairs to plot
source_sink_pair <- data.frame(source = c("psi","VPD","TA"),
                               sink = c("ET","ET","ET"))

# Make maps ===================
# Initialize a list to store continuous maps
g_continuous_ls <- list()
# Initialize a df to store Delta values
Delta_df <- c()
# Initialize a list to store TEtotal vs Delta scatter plots
g_TE_Delta_ls <- list()
# Initialize a list to store all maps (both continuous and categorical)
g_map_all_ls <- list()
# Initialize a list to store scatter and box plots
g_scatter_box_ls <- list()
# Initialize a list to store supplementary figures
g_scatter_box_SI_ls <- list()

# Colors for scatter and box plots
syc_colors <- RColorBrewer::brewer.pal(10,"Paired")[c(1,7,9)]
for(i in 1:nrow(source_sink_pair)){
  # Data processing -------------------
  source_name <- source_sink_pair$source[i]
  sink_name <- source_sink_pair$sink[i]
  
  # Get target synchrony metrics in the two directions
  syc_df1 <- read_syc_df(source_name,sink_name,Syc_metrics_path) %>%
    select(site_ID,syc1 = !!rlang::sym(res_name))
  syc_df2 <- read_syc_df(sink_name,source_name,Syc_metrics_path) %>%
    select(site_ID,syc2 = !!rlang::sym(res_name))
  
  # Merge df in the two directions
  syc_df <- merge(syc_df1,syc_df2,by="site_ID") %>%
    # Calculate delta = syc1 - syc2
    mutate(Delta = syc1 - syc2,
           source_sink = paste0(source_name,"_",sink_name)) %>%
    # Calculate total TE (syc1 + syc2)
    mutate(TE_total = syc1 + syc2,
           logTE_total = log(TE_total)) %>%
    # Calculate normalized asymmetry index
    mutate(A = Delta/TE_total) %>%
    # Merge with predictor df
    merge(predictor_df %>%
            rename(site_ID = site_id),by="site_ID") %>%
    # Remove NA
    na.omit() %>%
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
    inner_join(Site_info %>%
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
    Soil_Group = factor(Soil_Group,levels=c("Fine","Moderately fine","Medium","Moderately coarse","Coarse"))) %>%
    # Filter out outlier
    filter(site_ID!="US-CS1" & site_ID!= "US-NC1")

  # Store Delta synchrony metrics to get distribution
  Delta_df <- rbind(Delta_df,syc_df)
  
  # Map of continuous Delta, sized by continuous TE_total -----------------------
  # Continuous color palette for Delta
  Delta_color <- rev(RColorBrewer::brewer.pal(n=11,name = "PiYG"))
  g_map_continuous_Delta <- syc_map(syc_df,varname = "Delta",colors=Delta_color,
                                    legend_title = Delta_title,g_title = paste0(source_name,"-",sink_name),
                                    color_limits = c(-2.5,2.5),end_marks = c("left","right"),
                                    base_limits = c(0,2),base_raster = AI_raster,
                                    size_var = "TE_total",size_title = TEtotal_title)+
    theme(legend.position = "bottom",
          legend.box = "vertical")
  
  g_continuous_ls[[i]] <- g_map_continuous_Delta
  g_map_all_ls[[length(g_map_all_ls)+1]] <- g_map_continuous_Delta
  
  # Scatter plot of TE total vs Delta --------------------
  # Calculate thresholds to define regimes
  # This is the threshold for total TE, below which is defined as "Weakly coupled"
  TE_thr <- quantile(syc_df$TE_total,0.25,na.rm=TRUE)
  # Symmetry band around Delta = 0
  Delta_thr <- quantile(abs(syc_df$Delta),0.25,na.rm=TRUE)
  # Define regimes
  syc_df <- syc_df %>%
    mutate(
      regime = case_when(
        TE_total < TE_thr ~ "Weakly coupled",
        abs(Delta) <= Delta_thr ~ "Bidirectional",
        Delta > Delta_thr ~ "Driver-dominant",
        Delta < -Delta_thr ~ "ET-dominant"
      )
    )
  Regime_color <- RColorBrewer::brewer.pal(4,"Set2")
  g_Delta_TE_total_scatter_base <- ggplot(syc_df,
                                     aes(x = Delta, y = TE_total, color = regime)) +
    geom_point(size = 2.5, alpha = 0.7) +
    # Threshold lines
    geom_vline(xintercept = c(-Delta_thr, 0, Delta_thr),
               linetype = c("dashed","solid","dashed")) +
    geom_hline(yintercept = TE_thr, linetype = "dashed") +
    scale_color_manual(values = c(
      "Weakly coupled" = "grey70",
      "Bidirectional" = Regime_color[3],
      "Driver-dominant" = Regime_color[2],
      "ET-dominant" = Regime_color[1]
    )) +
    labs(
      x = Delta_title,
      y = TEtotal_title,
      color = ""
    ) +
    my_theme
  # Add violin plots for distribution
  dens_right <- axis_canvas(g_Delta_TE_total_scatter_base,axis = "y",coord_flip = TRUE)+
    geom_density(data=syc_df,aes(x=TE_total),
                 outline.type = "full",
                 fill = "grey",
                 alpha=0.3)+
    coord_flip()
  
  dens_top <- axis_canvas(g_Delta_TE_total_scatter_base, axis = "x") +
    geom_density(data = syc_df,aes(x = Delta),
                 outline.type = "full",fill = "grey",alpha = 0.3)
  
  g_Delta_TE_total_scatter <- g_Delta_TE_total_scatter_base
  
  g_Delta_TE_total_scatter <- insert_yaxis_grob(g_Delta_TE_total_scatter,
                                                dens_right,position = "right")

  g_Delta_TE_total_scatter <- insert_xaxis_grob(g_Delta_TE_total_scatter,
                                                dens_top,position = "top")
  
  g_TE_Delta_ls[[i]] <- g_Delta_TE_total_scatter
  
  # Map of categorical Delta -------------------- 
  regime_colors <- c(
    "Weakly coupled" = "grey70",
    "Bidirectional" = Regime_color[3],
    "Driver-dominant" = Regime_color[2],
    "ET-dominant" = Regime_color[1]
  )
  g_map_discrete_Delta <- syc_map_disc(syc_df,regime_var = "regime",
                                       regime_colors = regime_colors,
                                       base_raster = AI_raster,legend_title = "",
                                       base_limits = c(0,2),g_title = "")
    
  g_map_all_ls[[length(g_map_all_ls)+1]] <- g_map_discrete_Delta
  
  
  # Box plots and scatter plots of Delta across variable groups ---------------
  # scatter plot of Delta vs AI -----
  g_scatter <- scatter_vars(syc_df,"AI_gridded","Delta","source_sink",
                            "Aridity Index",Delta_title,syc_colors[i])+
    theme(legend.position = "none")
  g_scatter_box_SI_ls[[length(g_scatter_box_SI_ls)+1]] <- g_scatter
  
  # Box plot of Delta across aridity classes ------
  g_box_AI <- plot_box_groups(syc_df,"Delta","AI_Class",fill_color = syc_colors[i],
                              x_title = Delta_title,y_title = "",h_just = 0)
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_AI
  
  # Box plot of Delta across climate ----------
  # Aggregate climate classes
  g_box_climate_agg <- plot_box_groups(syc_df,"Delta","Koppen_aggregate",fill_color = syc_colors[i],
                                       x_title = Delta_title,y_title = "",h_just = 0)
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_climate_agg
  
  # Koppen climate without aggregation
  g_box_climate <- plot_box_groups(syc_df,"Delta","Koppen_clim_class",fill_color = syc_colors[i],
                                   x_title = Delta_title,y_title = "",box_violin = "Box",h_just = 0)
  g_scatter_box_SI_ls[[length(g_scatter_box_SI_ls)+1]] <- g_box_climate
  
  # Box plot of Delta across IGBP ---------
  # Aggregate IGBP classes
  g_box_IGBP_agg <- plot_box_groups(syc_df,"Delta","IGBP_aggregate",fill_color = syc_colors[i],
                                    x_title = Delta_title,y_title = "")
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_IGBP_agg
  
  # IGBP without aggregation
  g_box_IGBP <- plot_box_groups(syc_df,"Delta","IGBP_veg",fill_color = syc_colors[i],box_violin = "Box",
                                x_title = Delta_title,y_title = "")
  g_scatter_box_SI_ls[[length(g_scatter_box_SI_ls)+1]] <- g_box_IGBP
  
  # Box plots of Delta across soil type ------------
  # Soil with aggregation
  g_box_soil_agg <- plot_box_groups(syc_df,"Delta","Soil_Group",fill_color = syc_colors[i],
                                    x_title = Delta_title,y_title = "")
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_soil_agg
  
  # Soil without aggregation
  g_box_soil <- plot_box_groups(syc_df,"Delta","Soil_Type",fill_color = syc_colors[i],box_violin = "Box",
                                x_title = Delta_title,y_title = "")
  g_scatter_box_SI_ls[[length(g_scatter_box_SI_ls)+1]] <- g_box_soil
  
}

# Combine continuous Delta maps
g_continuous_map <- plot_grid(plotlist = g_continuous_ls,ncol=2,align="hv",labels = "auto")
#print_g(g_continuous_map,paste0("Delta_maps_AI_",res_name),14,10)
print_g(g_continuous_map,paste0("Delta_maps_AI_",res_name,"_LAI_filtered"),14,10)

# Combine continuous maps, and categorical maps
g_regimes <- plot_grid(plotlist = g_map_all_ls,ncol=2,labels="auto")
#print_g(g_regimes,paste0("Delta_regimes_",res_name),14,15)
print_g(g_regimes,paste0("Delta_regimes_",res_name,"_LAI_filtered"),14,15)

# Combine TEtotal vs Delta scatter plots
g_TE_Delta <- plot_grid(plotlist = g_TE_Delta_ls,ncol = 1,align="hv")
#print_g(g_TE_Delta,paste0("TEtotal_Delta_scatter_",res_name),4,12)
print_g(g_TE_Delta,paste0("TEtotal_Delta_scatter_",res_name,"_LAI_filtered"),4,12)
  
# Combine boxplots (Main text)
g_scatter_box <- plot_grid(plotlist = g_scatter_box_ls[c(1,5,9,2,6,10,3,7,11,4,8,12)],
                           ncol=3,align="hv",labels="auto")
#print_g(g_scatter_box,"Delta_compare_cli_veg_soil",16,9)
print_g(g_scatter_box,"Delta_compare_cli_veg_soil_LAI_filtered",16,9)

# Combine all scatter and box plots (SI figure)
g_scatter_box_SI <- plot_grid(plotlist = g_scatter_box_SI_ls[c(1,5,9,2,6,10,3,7,11,4,8,12)],
                              ncol=3,align="hv",labels = "auto")
#print_g(g_scatter_box_SI,"Delta_compare_cli_veg_soil_SI",16,14)
print_g(g_scatter_box_SI,"Delta_compare_cli_veg_soil_SI_LAI_filtered",16,14)


# Compare synchrony metrics across pairs =============
# Wilcoxon tests
my_comparisons <- list(
  c("psi_ET", "VPD_ET"),
  c("psi_ET", "TA_ET"),
  c("VPD_ET", "TA_ET")
)

# Process syc_df_tmp to make sure each site has all three Delta values
syc_df_tmp_paired <- Delta_df %>%
  select(site_ID,source_sink,value = Delta) %>%
  # Convert to long data to force matching
  pivot_wider(names_from = source_sink, values_from = value) %>%
  select(site_ID, psi_ET, VPD_ET, TA_ET) %>%
  na.omit() %>%
  # Convert back to long data
  pivot_longer(cols = c(psi_ET,VPD_ET,TA_ET),
               names_to = "source_sink",
               values_to = "value") %>%
  mutate(source_sink = factor(source_sink, levels = c("psi_ET","VPD_ET","TA_ET")))

# Compare synchrony metrics (ET as endpoint) across different pairs
g_box_pair <- ggplot(data = syc_df_tmp_paired,aes(x=source_sink,y=value,fill=source_sink,color=source_sink))+
  geom_half_violin(alpha=0.5,color=NA)+
  geom_boxplot(width=0.1,color="black",outlier.color = NA)+
  geom_jitter(aes(x = as.numeric(source_sink)+0.2),
              position = position_jitter(width=0.1),
              alpha=0.7,size=1)+
  my_theme+
  labs(y = Delta_title,x="")+
  scale_fill_manual(values = pair_color[c(1,3,5)])+
  scale_color_manual(values = pair_color[c(1,3,5)])+
  scale_x_discrete(labels = c(
    psi_ET = expression(psi %->% ET),
    VPD_ET = expression(VPD %->% ET),
    TA_ET  = expression(T[air] %->% ET)
  ))+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    paired = FALSE,
    label = "p.format",
    tip.length = 0.01,
    size=3
  )+
  # Set the top of the figure to make some room
  scale_y_continuous(expand = expansion(mult = c(0.05,0.15)))

#print_g(g_box_pair,paste0("Delta_box_",res_name),4,3)
print_g(g_box_pair,paste0("Delta_box_",res_name,"_LAI_filtered"),4,3)


summary(syc_df_tmp_paired$value[syc_df_tmp_paired$source_sink == "psi_ET"])
summary(syc_df_tmp_paired$value[syc_df_tmp_paired$source_sink == "VPD_ET"])
summary(syc_df_tmp_paired$value[syc_df_tmp_paired$source_sink == "TA_ET"])

