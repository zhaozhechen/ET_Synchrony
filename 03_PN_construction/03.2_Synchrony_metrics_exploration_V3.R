# Author: Zhaozhe Chen
# Date: 2026.2.8

# This code is to make maps for Delta

# -------- Global ----------
library(dplyr)
library(tidyr)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Predictor df
predictor_df <- read.csv("00_Data/perdictor_df_updated.csv")

# Source plotting functions
source("05_Visualization/Plotting_functions.R")
source("03_PN_construction/Synchrony_metrics_functions_v2.R")

Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Source-sink_pair-comparisons_V2/"

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
# Initialize a list to store the figures
g_ls <- list()
# Initialize a df to store Delta values
Delta_df <- c()
for(i in 1:nrow(source_sink_pair)){
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
    # Merge with predictor df
    merge(predictor_df %>%
            rename(site_ID = site_id),by="site_ID") %>%
    # Remove NA
    na.omit()
  
  # Get legend title
  #legend_expr <- bquote(
  #  TE(.(as.name(source_name)) %->% .(as.name(sink_name))) -
  #    TE(.(as.name(sink_name)) %->% .(as.name(source_name)))
  #)
  
  # Make map of Delta value
  g_map <- syc_map(syc_df,"Delta",palette_name = "PiYG",
                   g_title="",
                   legend_title="Delta",
                   color_limits=c(-2.5,2.5))  
  g_ls[[i]] <- g_map
  
  # Store Delta synchrony metrics to get distribution
  Delta_df <- rbind(Delta_df,syc_df)
}

# Combine maps
g_map_all <- plot_grid(plotlist = g_ls,
                       nrow=2,align="v",labels = "auto")

# Output this figure
print_g(g_map_all,paste0("Delta_maps_",res_name),
        12,6)

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



