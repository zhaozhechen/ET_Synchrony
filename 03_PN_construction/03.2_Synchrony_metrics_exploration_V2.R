# Author: Zhaozhe Chen
# Update Date: 2026.2.6

# This code is to debug previous results about psi -> ET vs ET -> psi

# -------- Global ----------
library(dplyr)
library(tidyr)

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
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/"

season <- "GS"

# Determine which metric to process
arrayid <- 1
# Name of target synchrony metrics
res_varname <- c("daily_p_TE","daily_agg_TE","best_lag")[arrayid]
# Title of these metrics
res_title <- c("Peak TE (%)","Daily aggregated TE (%)","Lag (h)")[arrayid]

my_color <- RColorBrewer::brewer.pal(5,"Set2")[c(1,2)]
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
g_map_ls <- list()
for(i in 1:nrow(sc_pairs)){
  source_name <- sc_pairs$source[i]
  sink_name <- sc_pairs$sink[i]
  # Remove NA for the target df
  syc_df_sub <- syc_df %>%
    filter(source_sink == paste0(source_name,"_",sink_name)) %>%
    filter(!is.na(.data[[res_name]]))
  # Make map
  # Continuous color palette
  map_color <- wes_palette("Zissou1",100,type="continuous")
  g_map <- syc_map(syc_df_sub,res_name,colors = map_color,
                   legend_title = "test",g_title = "",color_limits = c(0,15),end_marks = "right")
  g_map_ls[[i]] <- g_map
}
# Combine maps
g_map_all <- plot_grid(plotlist = g_map_ls,nrow=2,labels="auto")



# Synchrony metrics vs AI

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
print_g(g_box,paste0("Source-sink_pair-comparisons_V2/Box_",res_name),6,6)


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





