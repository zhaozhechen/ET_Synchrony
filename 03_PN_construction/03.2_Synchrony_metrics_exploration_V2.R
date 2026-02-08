# Author: Zhaozhe Chen
# Update Date: 2026.2.7

# This code is to explore synchrony strength among pairs

# -------- Global ----------
library(dplyr)
library(tidyr)
library(ggpubr)
library(FSA) # For Dunn test

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Predictor df, with updated AI_gridded
predictor_df <- read.csv("00_Data/perdictor_df_updated.csv")

# Source-sink pairs to test
source_sink_pairs <- read.csv("00_Data/Source-sink_pairs.csv")

# Source plotting functions
source("05_Visualization/Plotting_functions.R")
source("03_PN_construction/Synchrony_metrics_functions_v2.R")

# Output path of figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Source-sink_pair-comparisons_V2/"

season <- "GS"

# Determine which syc metric to process
arrayid <- 1
# Name of target synchrony metrics
res_varname <- c("daily_p_TE","daily_agg_TE","best_lag")[arrayid]
# Title of these metrics
res_title <- c("Peak TE (%)","Daily aggregated TE (%)","Lag (h)")[arrayid]

# This is the color pallet for paired comparisons
pair_color <- RColorBrewer::brewer.pal(10,"Paired")[c(1,2,7,8,9,10)]

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
  mutate(source_sink = factor(source_sink,levels = c("psi_ET","ET_psi","VPD_ET","ET_VPD","TA_ET","ET_TA")))

# Maps of synchrony metrics across pairs ==============================
# A list of source-sink pairs to plot
sc_pairs <- data.frame(source = c("psi","VPD","TA"),
                       sink = c("ET","ET","ET"))
# Initialize a list to store maps
g_map_ls <- list()
# Initialize a df to store synchrony subset
syc_df_tmp <- c()
# Initiliaze a list to store scatter and box plots
g_scatter_box_ls <- list()

# Colors for scatter and box plots
syc_colors <- RColorBrewer::brewer.pal(10,"Paired")[c(1,7,9)]
for(i in 1:nrow(sc_pairs)){
  source_name <- sc_pairs$source[i]
  sink_name <- sc_pairs$sink[i]
  # Remove NA for the target df
  syc_df_sub <- syc_df %>%
    filter(source_sink == paste0(source_name,"_",sink_name)) %>%
    filter(!is.na(.data[[res_name]]))
  
  # Store this df for later pairwise tests
  syc_df_tmp <- rbind(syc_df_tmp,syc_df_sub)
  
  # Make a scatter plot of target synchrony metric vs AI ---------
  g_scatter <- scatter_vars(syc_df_sub,"AI_gridded",res_name,"source_sink","Aridity Index",res_title,syc_colors[i])+
    theme(legend.position = "none")
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_scatter

  
  # Make a box plot of target synchrony metric across climate ----------
  g_box_climate <- plot_box(syc_df_sub,res_name,"Koppen_clim_class",fill_color = syc_colors[i],res_title,"")
  # Statistical test, compare target synchrony metric across climate group
  k_result <- kruskal.test(
    syc_df_sub[[res_name]] ~ syc_df_sub$Koppen_clim_class
  )
  # Get p-value
  p_txt <- paste0("p = ",formatC(k_result$p.value, format = "e", digits = 2))
  # Add p-value to the plot
  g_box_climate <- g_box_climate +
    annotate("text",label = p_txt,x=Inf,y=-Inf,hjust=1.1,vjust=-0.8,size=5)
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_climate  
  
  # Dunn test
  dunn_result <- dunnTest(x=syc_df_sub[[res_name]],
           g = syc_df_sub$Koppen_clim_class,
           method = "bh")$res
  
  # Make a box plot of target synchrony metric across IGBP ---------
  g_box_IGBP <- plot_box(syc_df_sub,res_name,"IGBP_veg",fill_color = syc_colors[i],res_title,"")
  # Statistical test, compare target synchrony metric across climate group
  k_result <- kruskal.test(
    syc_df_sub[[res_name]] ~ syc_df_sub$IGBP_veg
  )
  # Get p-value
  p_txt <- paste0("p = ",formatC(k_result$p.value, format = "e", digits = 2))
  # Add p-value to the plot
  g_box_IGBP <- g_box_IGBP +
    annotate("text",label = p_txt,x=Inf,y=-Inf,hjust=1.1,vjust=-0.8,size=5)
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_IGBP
  
  # Dunn test
  dunn_result <- dunnTest(x=syc_df_sub[[res_name]],
                          g = syc_df_sub$IGBP_veg,
                          method = "bh")$res
  
  # Make map
  # Continuous color palette
  map_color <- wes_palette("Zissou1",100,type="continuous")
  g_map <- syc_map(syc_df_sub,res_name,colors = map_color,
                   legend_title = res_title,g_title = "",color_limits = c(0,15),end_marks = "right")
  g_map_ls[[i]] <- g_map
}
# Combine 3 maps for three pairs
g_map <- plot_grid(plotlist = g_map_ls,nrow=2,align="hv",labels="auto")
print_g(g_map,paste0("Syc_maps_",res_name),12,6)

# Combine all scatter and box plots
g_scatter_box <- plot_grid(plotlist = g_scatter_box_ls,ncol=3,align="hv",labels="auto")
print_g(g_scatter_box,paste0("Syc_compare_",res_name),12,9)

# Compare synchrony metrics across pairs (for ET as endpoint) ==========
# Wilcoxon tests
my_comparisons <- list(
  c("psi_ET", "VPD_ET"),
  c("psi_ET", "TA_ET"),
  c("VPD_ET", "TA_ET")
)

# Process syc_df_tmp to make sure each site has three synchrony metrics
syc_df_tmp_paired <- syc_df_tmp %>%
  select(site_ID,source_sink,value = .data[[res_name]]) %>%
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
  labs(y = res_title,x="")+
  scale_fill_manual(values = pair_color[c(1,3,5)])+
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

print_g(g_box_pair,paste0("Syc_box_",res_name),4,3)

# Synchrony metrics vs AI ==================

# Synchrony metrics across climate

# Synchrony metrics across PFT


# Distributions of synchrony metrics across pairs ================================
# Box plots of synchrony metrics across different pairs 
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
print_g(g_box,paste0("Box_",res_name),6,6)


if(FALSE){
# Comparing target syc metric between two source-sink pairs =================
# Response synchrony metric name
res_varname <- "best_lag"

# Determines which source-sink pair to test
arrayid <- 1

source_name1 <- source_sink_pairs$Source1[arrayid]
sink_name1 <- source_sink_pairs$Sink1[arrayid]
source_name2 <- source_sink_pairs$Source2[arrayid]
sink_name2 <- source_sink_pairs$Sink2[arrayid]

# Make merged syc metrics, only keep target response variable (e.g., daily_p_TE), and calculate delta_response
pair_syc_df_merged <- merge_pair_syc_df(source_name1,sink_name1,source_name2,sink_name2,res_varname)

# Get title for the entire plot
g_title <- paste("Comparison of",res_varname,"between",
                 source_name1,"→",sink_name1,"vs",
                 source_name2,"→",sink_name2,"\n",
                 "+ =",res_varname,"(",source_name1,",",sink_name1,") >",
                 res_varname,"(",source_name2,",",sink_name2,")")

# Convert delta_response to binary values
pair_syc_df_merged_binary <- pair_syc_df_merged %>%
  mutate(delta_response = if_else(delta_response > 0,"+","-")) %>%
  # Remove NA
  filter(!is.na(delta_response)) %>%
  # Only keep growing season
  filter(Season == season)

# Make binary map
map_GS <- syc_map_disc(pair_syc_df_merged_binary %>% filter(Season == "GS"),
                       "delta_response","Delta",
                       g_title = "Growing Season")

# Check Sites with positive delta values (psi-> ET dominant than ET->psi)
# Make a map of these negative points, add site ID labels
g <- ggplot()+
  geom_sf(data=CONUS,fill="grey",color="black",alpha=0.3)+
  geom_point(data = pair_syc_df_merged_binary %>%
               filter(delta_response == "+"),
             aes(x=longitude,y=latitude),
             size=5,shape = 21,fill=my_color[1],color="black")+
  geom_text_repel(data = pair_syc_df_merged_binary %>%
                    filter(delta_response == "+"),
                  aes(x=longitude,y=latitude,label = site_ID))+
  map_theme
print_g(g,"Source-sink_pair-comparisons/daily_p_TE_psi_ET_greater_than_ET_psi",5,4)
}





