# Author: Zhaozhe Chen
# Update Date: 2026.3.31

# This code is to explore synchrony strength among pairs

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

# Determine which syc metric to process
arrayid <- 1
# Name of target synchrony metrics
res_varname <- c("daily_p_TE","daily_agg_TE","best_lag")[arrayid]
# Title of these metrics
res_title <- c("Peak TE (%)","Daily aggregated TE (%)","Lag (h)")[arrayid]

# This is the color pallet for paired comparisons
pair_color <- RColorBrewer::brewer.pal(10,"Paired")[c(1,2,7,8,9,10)]

# Colors for scatter and box plots
syc_colors <- RColorBrewer::brewer.pal(10,"Paired")[c(1,7,9)]

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

# Maps of synchrony metrics across pairs ==============================
# A list of source-sink pairs to plot
sc_pairs <- data.frame(source = c("psi","VPD","TA"),
                       sink = c("ET","ET","ET"))
# Initialize a list to store maps
g_map_ls <- list()
# Initialize a df to store synchrony subset
syc_df_tmp <- c()
# Initialize a list to store scatter and box plots
g_scatter_box_ls <- list()
# Initialize a list to store scatter plots and box plots for SI
g_scatter_box_SI_ls <- list()

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
  g_scatter_box_SI_ls[[length(g_scatter_box_SI_ls)+1]] <- g_scatter

  # Make a box plot of target synchrony metric across aridity classes -------
  g_box_AI <- plot_box_groups(syc_df_sub,res_name,"AI_Class",fill_color = syc_colors[i],
                              x_title = res_title,y_title = "",h_just = 0)
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_AI
  
  # Make a box plot of target synchrony metric across climate ----------
  # Aggregate climate classes
  g_box_climate_agg <- plot_box_groups(syc_df_sub,res_name,"Koppen_aggregate",fill_color = syc_colors[i],
                                       x_title = res_title,y_title = "",h_just = 0)
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_climate_agg
  
  # Koppen climate without aggregation
  g_box_climate <- plot_box_groups(syc_df_sub,res_name,"Koppen_clim_class",fill_color = syc_colors[i],
                                   x_title = res_title,y_title = "",box_violin = "Box",h_just = 0)
  g_scatter_box_SI_ls[[length(g_scatter_box_SI_ls)+1]] <- g_box_climate
  
  # Make a box plot of target synchrony metric across IGBP ---------
  # Aggregate IGBP classes
  g_box_IGBP_agg <- plot_box_groups(syc_df_sub,res_name,"IGBP_aggregate",fill_color = syc_colors[i],
                                    x_title = res_title,y_title = "")
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_IGBP_agg
  
  # IGBP without aggregation
  g_box_IGBP <- plot_box_groups(syc_df_sub,res_name,"IGBP_veg",fill_color = syc_colors[i],box_violin = "Box",
                                x_title = res_title,y_title = "")
  g_scatter_box_SI_ls[[length(g_scatter_box_SI_ls)+1]] <- g_box_IGBP
  
  # Make a box plot of target synchrony metric across soil type ------------
  # Soil with aggregation
  g_box_soil_agg <- plot_box_groups(syc_df_sub,res_name,"Soil_Group",fill_color = syc_colors[i],
                                    x_title = res_title,y_title = "")
  g_scatter_box_ls[[length(g_scatter_box_ls)+1]] <- g_box_soil_agg
  
  # Soil without aggregation
  g_box_soil <- plot_box_groups(syc_df_sub,res_name,"Soil_Type",fill_color = syc_colors[i],box_violin = "Box",
                                x_title = res_title,y_title = "")
  g_scatter_box_SI_ls[[length(g_scatter_box_SI_ls)+1]] <- g_box_soil

  # Make map  ----------------------------
  # Continuous color palette for synchrony values
  map_color <- RColorBrewer::brewer.pal(n=11,name = "YlGn")
  # Color limit: Use 0,10 for daily peak TE; Use 0.24 for lag
  # end_marks: use "right" for daily peak TE; Use "none" for lag
  g_map <- syc_map(syc_df_sub,res_name,colors = map_color,
                   legend_title = res_title,g_title = "",
                   color_limits = c(0,10),end_marks = "right",base_raster = AI_raster)+
    theme(legend.position = "bottom")
  
  
  g_map_ls[[i]] <- g_map
}

# Combine 3 maps for three pairs
g_map <- plot_grid(plotlist = g_map_ls,nrow=2,align="hv",labels="auto")
print_g(g_map,paste0("Syc_maps_AI_",res_name),14,7)

# Combine all scatter and box plots (main text)
g_scatter_box <- plot_grid(plotlist = g_scatter_box_ls[c(1,5,9,2,6,10,3,7,11,4,8,12)],
                           ncol=3,align="hv",labels="auto")
print_g(g_scatter_box,paste0("Syc_compare_cli_veg_soil_",res_name),16,9)

# Combine all scatter and box plots (SI figure)
g_scatter_box_SI <- plot_grid(plotlist = g_scatter_box_SI_ls[c(1,5,9,2,6,10,3,7,11,4,8,12)],
                              ncol=3,align="hv",labels = "auto")
print_g(g_scatter_box_SI,paste0("Syc_compare_cli_veg_soil_SI_",res_name),16,14)




# Compare synchrony metrics across pairs (for ET as endpoint) ==========
# Wilcoxon tests
my_comparisons <- list(
  c("psi_ET", "VPD_ET"),
  c("psi_ET", "TA_ET"),
  c("VPD_ET", "TA_ET")
)

# Process syc_df_tmp to make sure each site has three synchrony metrics
syc_df_tmp_paired <- syc_df_tmp %>%
  select(site_ID,source_sink,Regime,value = .data[[res_name]]) %>%
  # Convert to long data to force matching
  pivot_wider(names_from = source_sink, values_from = value) %>%
  #select(site_ID, psi_ET, VPD_ET, TA_ET) %>%
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
  labs(y = res_title,x="")+
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

# Compare synchrony metrics (ET as endpoint) energy vs water-limited across different pairs
# Calculate p-value per facet
p_df <- syc_df_tmp_paired %>%
  filter(is.finite(value), !is.na(Regime), !is.na(source_sink)) %>%
  group_by(source_sink) %>%
  wilcox_test(value ~ Regime) %>%
  mutate(
    p.format = scales::pvalue(p, accuracy = 0.001),
    group1 = "Water-limited",
    group2 = "Energy-limited"
  ) %>%
  add_xy_position(x = "Regime")          # gives xmin/xmax and default y.position

g_box_WvsE <- ggplot(data = syc_df_tmp_paired,aes(x=Regime,y=value,fill=source_sink))+
  geom_boxplot(outlier.color = "grey")+
  my_theme+
  labs(y = res_title,x="")+
  scale_fill_manual(values = pair_color[c(1,3,5)])+
  facet_wrap(~source_sink)+
  stat_pvalue_manual(p_df,label = "p.format",
                     tip.length = 0.01,size=3)+
  # Set the top of the figure to make some room
  scale_y_continuous(expand = expansion(mult = c(0.05,0.15)))

# Check outliers:
test <- syc_df_tmp_paired %>%
  arrange(source_sink,Regime,desc(value))

# Check sample size:
syc_df_tmp_paired %>%
  count(source_sink,Regime)

# Combine the two boxplots
g_box <- plot_grid(g_box_pair,g_box_WvsE,align="hv")

print_g(g_box,paste0("Syc_box_v2_",res_name),8,3)

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
print_g(g_box,paste0("Box_",res_name),3,4)

#plot_pdf(syc_df,res_name,"source_sink",pair_color,1,"test")


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


syc_df %>%
  group_by(source_sink) %>%
  summarize(mean = mean(GS_best_lag,na.rm=TRUE),
            median = median(GS_best_lag,na.rm=TRUE))




