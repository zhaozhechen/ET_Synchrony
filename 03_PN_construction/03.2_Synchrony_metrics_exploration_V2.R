# Author: Zhaozhe Chen
# Date: 2025.11.25

# This code is to debug previous results about psi -> ET vs ET -> psi

# -------- Global ----------
library(dplyr)
library(tidyr)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Predictor df
predictor_df <- read.csv("00_Data/perdictor_df.csv")

# Source-sink pairs to test
source_sink_pairs <- read.csv("00_Data/Source-sink_pairs.csv")

# Source plotting functions
source("05_Visualization/Plotting_functions.R")
source("03_PN_construction/Synchrony_metrics_functions_v2.R")

# Output path of figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/"

# Response synchrony metric name
res_varname <- "daily_p_TE"

season <- "GS"

# Determines which source-sink pair to test
arrayid <- 1

my_color <- RColorBrewer::brewer.pal(5,"Set2")[c(1,2)]

# ----- Main -------
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
