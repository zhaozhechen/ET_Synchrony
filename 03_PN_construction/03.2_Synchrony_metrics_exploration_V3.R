# Author: Zhaozhe Chen
# Date: 2026.2.5

# This code is to explore and analyze synchrony metrics

# -------- Global ----------
library(dplyr)
library(tidyr)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Predictor df
predictor_df <- read.csv("00_Data/perdictor_df.csv")

# Source plotting functions
source("05_Visualization/Plotting_functions.R")
source("03_PN_construction/Synchrony_metrics_functions_v2.R")

Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Source-sink_pair-comparisons_V2/"

# Season
season <- "GS"

# Response variable list
res_var_ls <- c("daily_p_TE","daily_agg_TE","best_lag")

# ------- Main --------
arrayid <- 1
# Target synchrony metric name
res_varname <- paste0(season,"_",res_var_ls[arrayid])

# Source and sink pairs to plot
source_sink_pair <- data.frame(source = c("psi","VPD","TA"),
                               sink = c("ET","ET","ET"))

# Initialize a list to store the figures
g_ls <- list()
# Initialize a df to store Delta values
Delta_df <- c()
for(i in 1:nrow(source_sink_pair)){
  source_name <- source_sink_pair$source[i]
  sink_name <- source_sink_pair$sink[i]
  
  # Get target synchrony metrics in the two directions
  syc_df1 <- read.csv(paste0(Syc_metrics_path,"Syc_metrics_df_",source_name,"_",sink_name,".csv")) %>%
    select(site_ID,syc1 = !!rlang::sym(res_varname))
  
  syc_df2 <- read.csv(paste0(Syc_metrics_path,"Syc_metrics_df_",sink_name,"_",source_name,".csv")) %>%
    select(site_ID,syc2 = !!rlang::sym(res_varname))
  
  # Merge df in the two directions
  syc_df <- merge(syc_df1,syc_df2,by="site_ID") %>%
    # Calculate delta = syc1 - syc2
    mutate(Delta = syc1 - syc2,
           source_sink = paste0(source_name,"_",sink_name)) %>%
    # Merge with coordinates
    merge(Site_info %>%
            select(site_ID = site_id,latitude,longitude),by="site_ID") %>%
    # Remove NA
    na.omit()
  
  # Get legend title
  legend_expr <- bquote(
    TE(.(as.name(source_name)) %->% .(as.name(sink_name))) -
      TE(.(as.name(sink_name)) %->% .(as.name(source_name)))
  )
  
  # Make map of Delta value
  g_map <- syc_map(syc_df,"Delta",palette_name = "RdYlBu",
                   g_title="",
                   legend_title=legend_expr,
                   color_limits=c(-2.5,2.5))  
  g_ls[[i]] <- g_map
  
  # Store Delta synchrony metrics to get distribution
  Delta_df <- rbind(Delta_df,syc_df %>% select(Delta,source_sink))
}

# Make distribution plots of Delta
syc_3_colors <- brewer.pal(10,"Set3")[c(1,3,6)]
g_hist <- ggplot(data=Delta_df,aes(x = Delta,color=source_sink))+
  geom_density(fill=NA,linewidth=1)+
  my_theme+
  theme(legend.position = c(0.8,0.8))+
  labs(color="")+
  scale_color_manual(values = c("psi_ET" = syc_3_colors[2],
                                "VPD_ET" = syc_3_colors[1],
                                "TA_ET" = syc_3_colors[3]))

# Combine maps
g_map_all <- plot_grid(plotlist = g_ls,
                       ncol=1,align="v",labels = "auto")
# Output this figure
print_g(g_map_all,paste0("Delta_maps_",res_varname),
        10,12)







