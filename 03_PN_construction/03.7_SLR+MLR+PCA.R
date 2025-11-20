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

# arrayid determines the source and sink variable 
arrayid <- 1
# season includes "GS","NGS","FT". Determines which season to consider
season <- "GS"
# target_syc_metric_name determines which synchrony metric to consider
#target_syc_metric_name <- "daily_agg_TE"
target_syc_metric_name <- "daily_p_TE"
# The corresponding label for this target syc metric
#target_syc_metric_label <- bquote(DailyTE~"("~psi~"→"~ET~")")
target_syc_metric_label <- bquote(TE[max])

# Predictors to include
pre_var_ls <- c("AVG_SILT","AVG_CLAY","elevation","porosity","CH_GLAD","AI","RD")

# Labels for all variables
my_label <- c(target_syc_metric_label,"Silt (%)","Clay (%)","Elevation",
              "Porosity","Canopy height","Aridity index","Rooting depth")

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
CM <- Hmisc::rcorr(as.matrix(syc_df[,-1]),type = "spearman")
g_CM <- plot_CM(CM)
print_g(g_CM,paste0("CM_",season,"_",target_syc_metric_name,"_",source_name,"_",sink_name),8,8)




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






