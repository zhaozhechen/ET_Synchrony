# Author: Zhaozhe Chen
# Date: 2026.3.20

# This code is to make maps for Delta

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
    na.omit()

  # Store Delta synchrony metrics to get distribution
  Delta_df <- rbind(Delta_df,syc_df)
  
  # Map of continuous Delta, sized by continuous TE_total -----------------------
  # Continuous color palette for Delta
  Delta_color <- rev(RColorBrewer::brewer.pal(n=11,name = "PiYG"))
  g_map_continuous_Delta <- syc_map(syc_df,varname = "Delta",colors=Delta_color,
                                    legend_title = "Directional asymmetry (ΔTE, %)",g_title = paste0(source_name,"-",sink_name),
                                    color_limits = c(-2.5,2.5),end_marks = c("left","right"),
                                    base_limits = c(0,2),base_raster = AI_raster,
                                    size_var = "TE_total",size_title = "Bidirectional coupling strength (TEtotal, %)")+
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
      x = "Directional asymmetry (ΔTE, %)",
      y = "Bidirectional coupling strength (TEtotal, %)",
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
  
}

# Combine continuous Delta maps
g_continuous_map <- plot_grid(plotlist = g_continuous_ls,ncol=2,align="hv",labels = "auto")
print_g(g_continuous_map,paste0("Delta_maps_AI_",res_name),14,10)
# Combine continuous maps, and categorical maps
g_regimes <- plot_grid(plotlist = g_map_all_ls,ncol=2,labels="auto")
print_g(g_regimes,paste0("Delta_regimes_",res_name),14,15)
# Combine TEtotal vs Delta scatter plots
g_TE_Delta <- plot_grid(plotlist = g_TE_Delta_ls,ncol = 1,align="hv")
print_g(g_TE_Delta,paste0("TEtotal_Delta_scatter_",res_name),4,12)
  
  
  

  
  # Make a scatter plot of Delta vs AI ----------------
  g_scatter <- scatter_vars(syc_df,"AI_gridded","Delta","source_sink","Aridity Index",
                            "Directional difference in TE(%)",syc_colors[i])+
    theme(legend.position = "none")
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_scatter
  
  # Make a box plot of Delta across climate -----------
  g_box_climate <- plot_box(syc_df,"Delta","Koppen_clim_class",fill_color=syc_colors[i],
                            "Directional difference in TE(%)","")
  # Statistical test, compare target synchrony metric across climate group
  k_result <- kruskal.test(
    syc_df$Delta ~ syc_df$Koppen_clim_class
  )
  # Get p-value
  p_txt <- paste0("p = ",formatC(k_result$p.value, format = "e", digits = 2))
  # Add p-value to the plot
  g_box_climate <- g_box_climate +
    annotate("text",label = p_txt,x=Inf,y=-Inf,hjust=1.1,vjust=-0.8,size=5)
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_climate  
  
  # Make a box plot of Delta across IGBP ---------
  g_box_IGBP <- plot_box(syc_df,"Delta","IGBP_veg",fill_color = syc_colors[i],
                         "Directional difference in TE(%)","")
  # Statistical test, compare target synchrony metric across climate group
  k_result <- kruskal.test(
    syc_df$Delta ~ syc_df$IGBP_veg
  )
  # Get p-value
  p_txt <- paste0("p = ",formatC(k_result$p.value, format = "e", digits = 2))
  # Add p-value to the plot
  g_box_IGBP <- g_box_IGBP +
    annotate("text",label = p_txt,x=Inf,y=-Inf,hjust=1.1,vjust=-0.8,size=5)
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_IGBP
  

}

# Combine maps
g_map_all <- plot_grid(plotlist = g_ls,
                       nrow=2,align="v",labels = "auto")

# Output this figure
print_g(g_map_all,paste0("Delta_maps_",res_name),
        12,6)

# Combine all scatter and box plots
g_scatter_box <- plot_grid(plotlist = g_scatter_box_ls,ncol=3,align="hv",labels="auto")
print_g(g_scatter_box,paste0("Delta_compare_",res_name),12,9)

# Compare synchrony metrics across pairs =============
# Wilcoxon tests
my_comparisons <- list(
  c("psi_ET", "VPD_ET"),
  c("psi_ET", "TA_ET"),
  c("VPD_ET", "TA_ET")
)

# Process syc_df_tmp to make sure each site has three synchrony metrics
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

g_box_pair <- ggplot(data = syc_df_tmp_paired,aes(x=source_sink,y=value,fill=source_sink))+
  geom_boxplot(outlier.color = "grey")+
  my_theme+
  labs(y = paste0("Delta ",res_title),x="")+
  scale_fill_manual(values = pair_color[c(1,3,5)])+
  scale_x_discrete(labels = c(
    psi_ET = expression(psi~","~ET),
    VPD_ET = expression(VPD~","~ET),
    TA_ET  = expression(T[air]~","~ET)
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

print_g(g_box_pair,paste0("Delta_box_",res_name),4,3)



