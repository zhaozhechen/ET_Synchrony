# Author: Zhaozhe Chen
# Date: 2025.11.19

# This code is to explore how CZ structures affect synchrony metrics
# Using approaches including Correlation matrix (CM), SLR, MLR, and PCA

# -------- Global --------
library(dplyr)
library(factoextra)
library(FactoMineR)
library(ggrepel)
library(vegan)
library(corrplot)
library(Hmisc) # this is for correlation matrix
library(reshape2)
library(RColorBrewer)

# Data path ==============================
# Path to synchrony metrics df 
syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"
# Read predictor df (for AI)
predictor_df <- read.csv("00_Data/perdictor_df_updated.csv")
# Output path for figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_vs_predictors_CM+SLR+MLR+PCA/"

# Source Functions ===========
source("05_Visualization/Plotting_functions.R")

# Global variables to determine what metrics to process ============
# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")
# All variable combinations
var_comb <- expand.grid(from = var_ls,to = var_ls) %>%
  filter(from != to)

# List of synchrony metric names
target_syc_metric_name_ls <- c("daily_agg_TE","daily_p_TE")
# List of their name labels
target_syc_metric_label_ls <- list(bquote(DailyTE~"("~psi~"→"~ET~")"),bquote(TE[max]))

# arrayid determines the source and sink variable 
arrayid <- 1
# season includes "GS","NGS","FT". Determines which season to consider
season <- "GS"
# Determine which syc metric to process
syc_id <- 1
target_syc_metric_name <- target_syc_metric_name_ls[syc_id]
target_syc_metric_label <- target_syc_metric_label_ls[[syc_id]]

# Predictors to include
pre_var_ls <- c("AVG_SILT","AVG_CLAY","elevation","porosity","CH_GLAD","AI","RD",
                "avg_air_temp_degC","avg_precip_mm","Koppen_clim_class","IGBP_veg")

# Labels for all predictors
pre_var_label_ls <- list("Silt (%)","Clay (%)","Elevation","Porosity","Canopy height","Aridity index","Rooting depth",
                      bquote("Long-term"~T[air]),"Long-term P","Koppen Climate Class","PFT")

# ------ Main ----------
# Step 1. Process synchrony metrics and match with predictors ============
# Extract target synchrony metric values
source_name <- var_comb$from[arrayid]
sink_name <- var_comb$to[arrayid]
# Read in synchrony metrics df between these two pairs
syc_df <- read.csv(paste0(syc_metrics_path,"Syc_metrics_df_",source_name,"_",sink_name,".csv"))
# Extract only the target synchrony metric
syc_df <- syc_df %>%
  select(site_ID,paste0(season,"_",target_syc_metric_name)) %>%
  rename(site_id = site_ID)
# Get required predictors and merge with response (synchrony metric)
syc_df <- syc_df %>%
  left_join(predictor_df %>%
              select(c("site_id",pre_var_ls)),
            by="site_id")
  
# Step 2. Correlation matrix ==============
# Remove categorical data
CM <- Hmisc::rcorr(as.matrix(syc_df[,-c(1,12,13)]),type = "spearman")
# Labels for the variables
my_label <- c(target_syc_metric_label,pre_var_label_ls[-c(10,11)])
g_CM <- plot_CM(CM,my_label)
print_g(g_CM,paste0("CM_",season,"_",target_syc_metric_name,"_",source_name,"_",sink_name),10,10)

# Remove highly correlated variables from the predictors (r > 0.65)
pre_var_ls <- pre_var_ls[-c(1,4,9)]
pre_var_label_ls <- pre_var_label_ls[-c(1,4,9)]

# Step 3. SLR of syc metric vs predictors ================
# predictor names, excluding categorical data for SLR
pre_var_ls_SLR <- pre_var_ls[-c(7,8)]
# Corresponding name labels
pre_var_label_ls_SLR <- pre_var_label_ls[-c(7,8)]

# Initialize a list to store figures
g_scatter_ls <- list()
for(i in 1:length(pre_var_ls_SLR)){
  # Get predictor name
  pred_name <- pre_var_ls_SLR[i]
  # Its label name
  pred_label <- pre_var_label_ls_SLR[[i]]
  
  # Make scatter plots
  g_scatter <- scatter_vars(syc_df,pred_name,paste0(season,"_",target_syc_metric_name),
                            group_name = "IGBP_veg",xtitle = pred_label,ytitle = target_syc_metric_label,
                            my_color = brewer.pal(11,"Set3"))+
    theme(legend.position = "top")
  g_scatter_ls[[i]] <- g_scatter
}
g_scatter <- plot_grid(plotlist = g_scatter_ls,nrow=2,align = "hv")
print_g(g_scatter,paste0("Scatter_",season,"_",target_syc_metric_name,"_",source_name,"_",sink_name),10,8)

# Step 4. MLR of syc metric vs predictors ================


#MLR_fit <- lm(Y_scaled~X_scaled)





# Run PCA: predictors = active, synchrony = supplementary quantitative variables
#res.pca <- PCA(
#  merged_df[, c(predictor_vars, synchro_vars)],
#  scale.unit = TRUE,
#  ncp = 5,
#  graph = FALSE,
#  quanti.sup = (length(predictor_vars) + 1):(length(predictor_vars) + length(synchro_vars))
#)



# ---- Visualization ----

#fviz_pca_var(res.pca, repel = TRUE, col.var = "black") +
#  ggtitle("Predictor PCA with responses as supplementary variables")






