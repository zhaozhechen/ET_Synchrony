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

Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Source-sink_pair-comparisons_V3/"

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
# Initialize a list to store continuous figures
g_continuous_ls <- list()
# Initialize a df to store Delta values
Delta_df <- c()
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
                                    legend_title = "Delta",g_title = paste0(source_name,"-",sink_name),
                                    color_limits = c(-2.5,2.5),end_marks = c("left","right"),
                                    base_raster = AI_raster,
                                    size_var = "TE_total",size_title = "Bidirectional coupling strength")+
    theme(legend.position = "bottom")
  g_continuous_ls[[i]] <- g_map_continuous_Delta

  # Scatter plot of TE total vs Delta ------
  
  
  # Map of categorical Delta ---------- 
  
  
  
  
  # Make categorical map of Delta values ---------------
  # Calculate quantiles to divide categories
  q25 <- quantile(syc_df$Delta,probs = 0.25)
  q75 <- quantile(syc_df$Delta,probs = 0.75)
  syc_df_categorical <- syc_df %>%
    mutate(Type = case_when(
      Delta >= q75 ~ "Driver-dominant",
      Delta <= q25 ~ "ET-dominant",
      TRUE ~ "Near-symmetric"))
  g_map_categorical <- ggplot()+
    geom_sf(data = CONUS, fill = "grey", color = "black", alpha = 0.3) +
    geom_point(data = syc_df_categorical,aes(x=longitude.x,y=latitude.x,fill=Type),
               size=5,alpha=0.8,shape=21,color="black")+
    map_theme+
    scale_fill_manual(values = syc_colors)+
    labs(fill="")
  
  
  
  
  

  
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



