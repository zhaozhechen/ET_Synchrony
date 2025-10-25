# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.10.24

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

# Output path for figures
#Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_vs_predictors/"
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Source-sink_pair-comparisons/"

# Colors for the three seasons
season_color <- brewer.pal(3,"Set2")
# Palette for making maps
palette_name <- "YlOrRd"

# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")
# All seasons
seasons <- c("FT","GS","NGS")

# All variable combinations
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# Response variable list
res_var_ls <- c("daily_p_TE","daily_agg_TE","best_lag","mem1","mem2","mem3","mem4","mem5")
# Predictor variable list
pre_var_ls <- c("AI","CH","RD","TSand","elevation","porosity")

# Determine whether to make plots
If_plot <- FALSE

# --------- Main ---------
source_name1 <- "psi"
sink_name1 <- 'ET'

source_name2 <- "ET"
sink_name2 <- "psi"

# target response variable name
res_varname <- "daily_p_TE"

# Compare syc metrics across source-sink pairs =============
# Making plots for delta of the target variable (e.g., daily_p_TE) between two source-sink pairs

# Get title for the entire plot
g_title <- paste("Comparison of",res_varname,"between",
                 source_name1,"→",sink_name1,"vs",
                 source_name2,"→",sink_name2,"\n",
                 "Delta =",res_varname,"(",source_name1,",",sink_name1,") -",
                 res_varname,"(",source_name2,",",sink_name2,")")

# Make merged syc metrics, only keep target response variable (e.g., daily_p_TE), and calculate delta_response
pair_syc_df_merged <- merge_pair_syc_df(source_name1,sink_name1,source_name2,sink_name2,res_varname)

# Make plots for continuous values -----------------
# Maps of Delta (difference in response variable between the two source-sink pairs) across seasons
map_FT <- syc_map(pair_syc_df_merged %>% filter(Season == "FT"),
                     "delta_response",palette_name = "RdYlBu","Delta",
                     g_title = "Full Time Series",color_limits = c(-2,2))
map_GS <- syc_map(pair_syc_df_merged %>% filter(Season == "GS"),
                     "delta_response",palette_name = "RdYlBu","Delta",
                     g_title = "Growing Season",color_limits = c(-2,2))
map_NGS <- syc_map(pair_syc_df_merged %>% filter(Season == "NGS"),
                     "delta_response",palette_name = "RdYlBu","Delta",
                     g_title = "Non-Growing Season",color_limits = c(-2,2))

# Combine these maps
g_maps <- plot_grid(map_FT,map_GS,map_NGS,align = "hv",labels = "auto",nrow=1)

# Initialize a list to store pdf and scatter plots
g_ls <- list()
# Get distributions of differences (delta) of response variable between the two source-sink pairs, across seasons
g_ls[[1]] <- plot_pdf(pair_syc_df_merged,"delta_response","Season",
                      season_color,x_intercept = 0, x_title = "Delta")

# Make scatter plots across predictors
for(pre_varname in pre_var_ls){
  g_scatter <- syc_scatter_long(pair_syc_df_merged,"delta_response",pre_varname,y_title = "Delta")
  g_ls[[length(g_ls)+1]] <- g_scatter
}

# Combine pdf and scatter plots
g_scatters <- plot_grid(plotlist = g_ls,ncol=4,align = "hv",labels = letters[4:10])
# Combine with maps
g_all <- plot_grid(g_maps,g_scatters,align = "hv",
                   nrow=2,
                   rel_heights = c(1,2))
# Add title to the entire plot
g_all_titled <- plot_grid(
  ggdraw() + 
    draw_label(g_title, 
               fontface = "bold", 
               size = 14, 
               x = 0.5, y = 0.5, 
               hjust = 0.5, vjust = 0.5),
  g_all,
  ncol = 1,
  rel_heights = c(0.08, 1) 
)

# Output this figure
print_g(g_all_titled,paste0("Comparison_between_",source_name1,"-",sink_name1,"_vs_",source_name2,"-",sink_name2,"_continuous"),
        12,9)

# Make plots for binary values-> Focus on direction (sign) -------------------------
# Convert delta_response to binary values
pair_syc_df_merged_binary <- pair_syc_df_merged %>%
  mutate(delta_response = if_else(delta_response > 0,"+","-")) %>%
  # Remove NA
  filter(!is.na(delta_response))
# This is for positive and negative colors
my_color <- brewer.pal(11,"RdBu")[c(10,2)]
# Maps of Delta (difference in response variable between the two source-sink pairs) across seasons
map_FT <- syc_map_disc(pair_syc_df_merged_binary %>% filter(Season == "FT"),
                  "delta_response","Delta",
                  g_title = "Full Time Series")
map_GS <- syc_map_disc(pair_syc_df_merged_binary %>% filter(Season == "GS"),
                       "delta_response","Delta",
                       g_title = "Growing Season")
map_NGS <- syc_map_disc(pair_syc_df_merged_binary %>% filter(Season == "NGS"),
                       "delta_response","Delta",
                       g_title = "Non-Growing Season")

# Combine these maps
g_maps <- plot_grid(map_FT,map_GS,map_NGS,align = "hv",labels = "auto",nrow=1)

# Bar plot showing the percentage of positive vs negative during three seasons
g_prop_bar <- pair_syc_df_merged_binary %>%
  count(Season, delta_response) %>%
  group_by(Season) %>%
  mutate(Prop = n / sum(n)) %>%
  ggplot(aes(x = Season, y = Prop, fill = delta_response)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = my_color, name = "Direction") +
  labs(y = "Proportion", x = "Season") +
  my_theme +
  theme(legend.position = "top")

# Proportion of positive (+) values across predictor bins
for(pre_varname in pre_var_ls){
  g_prop_line <- prop_sign_plot(pair_syc_df_merged_binary,pre_varname,season_color)
}









pair_syc_df_merged_binary %>%
  mutate(AI_bin = cut(AI, breaks = quantile(AI, probs = seq(0, 1, 0.25), na.rm = TRUE))) %>%
  group_by(AI_bin, Season) %>%
  summarise(Prop_pos = mean(delta_response == "+", na.rm = TRUE)) %>%
  ggplot(aes(x = AI_bin, y = Prop_pos, group = Season, color = Season)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(x = "Aridity Index (quartiles)", y = "Proportion of + (ψ→ET)") +
  my_theme

ggplot(pair_syc_df_merged_binary, aes(x = delta_response, y = AI, fill = delta_response)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~ Season) +
  scale_fill_manual(values = my_color, name = "Direction") +
  labs(x = "Direction", y = "Aridity Index") +
  my_theme

glm_dir <- glm(I(delta_response == "+") ~ AI + CH + RD + TSand + elevation + porosity, 
               data = pair_syc_df_merged_binary, 
               family = binomial)

pair_syc_df_merged_binary %>%
  select(Site, Season, delta_response) %>%
  spread(Season, delta_response) %>%
  count(FT, GS, NGS)

# Scatter plots of syc metrics vs predictors for each of the 12 source-sink pairs ================
if(If_plot){
for(pre_id in 1:length(pre_var_ls)){
  # Variable name for predictor
  pre_varname <- pre_var_ls[pre_id]
  # Loop over each pair of source and sink variables
  # Initialize a list to store the figures
  g_ls <- list()
  for(pair_id in 1:12){
    # Determine which variable pair to process
    source_name <- var_comb$from[pair_id]
    sink_name <- var_comb$to[pair_id]
    # Read in synchorny metrics df
    Syc_df <- read.csv(paste0(Syc_metrics_path,"Syc_metrics_df_",source_name,"_",sink_name,".csv")) %>%
      dplyr::select(- X)
    
    # Merge with predictors
    Syc_df_merge <- Syc_df %>%
      # Merge with predictors
      left_join(predictor_df %>% dplyr::select(-X),by="site_ID") %>%
      # Merge with site info
      left_join(Site_info %>% dplyr::select(-X),by=c("site_ID" = "site_id"))
    
    # Loop over the response variables
    for(res_id in 1:length(res_var_ls)){
      res_varname <- res_var_ls[res_id]
    
      # Make a figure
      g_scatter <- syc_scatter(Syc_df_merge,res_varname,pre_varname)+
        ggtitle(paste(source_name,"→",sink_name))
      # Store this figure in the figure list
      g_ls[[length(g_ls)+1]] <- g_scatter
    }
  }
  
  # Combine all figures for this predictor into 1 plot
  g_all <- plot_grid(plotlist = g_ls,ncol = length(res_var_ls),
                     align = "hv")
  # Output this figure
  print_g(g_all,paste0("Syc_metrics_vs_",pre_varname),
          24,36)
  message(pre_id,"out of",length(pre_var_ls))
}
}

# Previous plotting codes =======================================================
if(FALSE){
  # Compare 5 memory options ========
  # Loop each season and variable pair
  for(i in 1:nrow(var_comb)){
    source_name <- var_comb$from[i]
    sink_name <- var_comb$to[i]
    # Get memory columns names
    mem_names <- paste0(c("mem1","mem2","mem3","mem4","mem5"),"_",source_name,"_to_",sink_name)
    # Initiliate a list to store three seasons
    g_seasons <- list()
    for(season in seasons){
      # Get a subset of these memory columns for this season
      df_tmp <- Syc_metrics_df %>%
        filter(GS == season) %>%
        select(Site_ID,all_of(mem_names))
      # Rename the columns
      names(df_tmp) <- c("Site_ID","Mem1","Mem2","Mem3","Mem4","Mem5")
      # Pivot to long data
      df_tmp_long <- df_tmp %>%
        pivot_longer(cols = starts_with("Mem"),
                     names_to = "Metric",
                     values_to = "value")
      # Create scatter plot matrix
      g_matrix <- ggpairs(df_tmp[,c(2:6)])+
        ggtitle(paste0("Memory comparison (",source_name," → ",sink_name,"; ",season,")"))+
        my_theme
      g_seasons[[season]] <- grob_to_ggplot(g_matrix)
    }
    
    
    # Output plots for this variable pair
    g <- wrap_plots(g_seasons,nrow=1)
    print_g(g,paste0("Mem5_comparison/Mem4_",source_name,"_to_",sink_name),
            15,5)
    
    message(i)
  }
  
  # Preprocessing of synchrony df -------------------
  Syc_metrics_df <- Syc_metrics_df  %>%
    select(-X) %>%
    # Join by Site_info
    left_join(Site_info %>% rename(Site_ID = site_id,Soil_type = Description),
              by = "Site_ID") %>%
    # Join by AI
    left_join(Site_info_AI %>% select(Site_ID,AI),by="Site_ID") %>%
    mutate(Soil_type = as.factor(Soil_type),
           IGBP_veg = as.factor(IGBP_veg),
           GS = as.factor(GS),
           Koppen_clim_class = as.factor(Koppen_clim_class)) %>%
    # Classify AI into five gradients
    mutate(
      AI_level = case_when(
        AI < 0.05 ~ "Hyperarid",
        AI >= 0.05 & AI < 0.2 ~ "Arid",
        AI >= 0.2 & AI < 0.5 ~ "Semiarid",
        AI >= 0.5 & AI < 0.65 ~ "Semihumid",
        AI >= 0.65 ~ "Humid"
      ),
      AI_level = factor(AI_level,
                        levels=c("Hyperarid","Arid","Semiarid","Semihumid","Humid"),
                        ordered = TRUE)
    )
  
  # Make plots for target variable ---------------------
  
  # List of syc metrics to plot
  varname_ls <- names(Syc_metrics_df)[2:61]
  
  # Loop over each variable
  for(varname in varname_ls){
    plot_syc_all(Syc_metrics_df,varname,palette_name)  
  }  
}













