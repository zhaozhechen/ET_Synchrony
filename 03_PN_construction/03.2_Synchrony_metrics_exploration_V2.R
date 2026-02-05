# Author: Zhaozhe Chen
# Update Date: 2025.12.11

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

my_color <- RColorBrewer::brewer.pal(5,"Set2")[c(1,2)]

# ----- Main -------
# Distributions of target syc metric =================================
# Name of target variable
res_varname <- "daily_agg_TE"
res_name <- paste0(season,"_",res_varname)

source_name1 <- "psi"
sink_name1 <- "ET"
source_name2 <- "VPD"
sink_name2 <- "ET"
source_name3 <- "ET"
sink_name3 <- "psi"

syc_df1 <- read.csv(paste0(Syc_metrics_path,"Syc_metrics_df_",source_name1,"_",sink_name1,".csv")) %>%
  mutate(source_sink = paste0(source_name1,"_",sink_name1))
syc_df2 <- read.csv(paste0(Syc_metrics_path,"Syc_metrics_df_",source_name2,"_",sink_name2,".csv")) %>%
  mutate(source_sink = paste0(source_name2,"_",sink_name2))
syc_df3 <- read.csv(paste0(Syc_metrics_path,"Syc_metrics_df_",source_name3,"_",sink_name3,".csv")) %>%
  mutate(source_sink = paste0(source_name3,"_",sink_name3))
  
syc_df <- rbind(syc_df1,syc_df2,syc_df3) %>%
  merge(predictor_df %>%
          rename(site_ID = site_id),by="site_ID")

# Get the mean and median of target syc metric
syc_mean1 <- mean(syc_df1[[res_name]],na.rm=TRUE)
syc_median1 <- median(syc_df1[[res_name]],na.rm=TRUE)

syc_mean2 <- mean(syc_df2[[res_name]],na.rm=TRUE)
syc_median2 <- median(syc_df2[[res_name]],na.rm=TRUE)

syc_mean3 <- mean(syc_df3[[res_name]],na.rm=TRUE)
syc_median3 <- median(syc_df3[[res_name]],na.rm=TRUE)

syc_3_colors <- brewer.pal(10,"Set3")[c(1,3,6)]
# Get distribution of target variable
g_hist <- ggplot(data=syc_df,aes(x=.data[[res_name]],color=source_sink))+
  geom_density(fill=NA,linewidth=1)+
  #geom_histogram()+
  my_theme+
  theme(legend.position = "right")+
  labs(color="",x=res_varname)+
  scale_color_manual(values = c("ET_psi" = syc_3_colors[1],
                                "psi_ET" = syc_3_colors[2],
                                "VPD_ET" = syc_3_colors[3]))+
  geom_vline(data = data.frame(source_sink = c(paste0(source_name1,"_",sink_name1),
                                               paste0(source_name2,"_",sink_name2),
                                               paste0(source_name3,"_",sink_name3)),
                               medians = c(syc_median1,syc_median2,syc_median3),
                               means = c(syc_mean1,syc_mean2,syc_mean3)),
             aes(xintercept = means,color=source_sink),
             linetype = "dashed",
             linewidth=1)
print_g(g_hist,paste0("Source-sink_pair-comparisons/",res_name,"_distributions"),5.5,4)




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






