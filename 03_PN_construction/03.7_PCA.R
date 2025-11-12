# Author: Zhaozhe Chen
# Date: 2025.11.11

# This code is to conduct PCA for synchrony metrics

# -------- Global --------
library(dplyr)
library(factoextra)
library(FactoMineR)
library(ggrepel)
library(vegan)
library(corrplot)

# Path to synchrony metrics df 
syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"
# Read predictor df (for AI)
predictor_df <- read.csv("00_Data/perdictor_df_updated.csv")

# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")
# All variable combinations
var_comb <- expand.grid(from = var_ls,to = var_ls) %>%
  filter(from != to)

season <- "GS"
target_syc_metric_name <- "daily_agg_TE"

# Predictors to include
pre_var_ls <- c("AVG_SILT","AVG_CLAY","elevation","porosity","CH_GLAD","AI","RD")

# ------ Main ----------
# Make combined syc metrics df (one column per pair)
# Initialize a list to store all df
syc_metrics_ls <- list()
for(arrayid in 1:nrow(var_comb)){
  source_name <- var_comb$from[arrayid]
  sink_name <- var_comb$to[arrayid]
  # Read in synchrony metrics df between these two pairs
  syc_df <- read.csv(paste0(syc_metrics_path,"Syc_metrics_df_",source_name,"_",sink_name,".csv"))
  # Get Site_ID once
  if(arrayid == 1){Site_ID_vec <- syc_df$site_ID}
  # Extract the target synchrony metric value
  syc_metrics <- syc_df[paste0(season,"_",target_syc_metric_name)]
  syc_metrics_ls[[arrayid]] <- syc_metrics
}

syc_metrics_df <- do.call(cbind,syc_metrics_ls)
names(syc_metrics_df) <- paste0(var_comb$from,"_to_",var_comb$to)
syc_metrics_df <- cbind(Site_ID_vec,syc_metrics_df)
names(syc_metrics_df)[1] <- "site_id"

# Get required predictor variables and merge with response (synchrony metrics ) --------------
pre_df <- predictor_df %>%
  select(c("site_id",pre_var_ls))

merged_df <- merge(pre_df,syc_metrics_df,by="site_id")

# Identify variable sets
predictor_vars <- pre_var_ls
synchro_vars <- grep("_to_", names(merged_df), value = TRUE)

# Run PCA: predictors = active, synchrony = supplementary quantitative variables
res.pca <- PCA(
  merged_df[, c(predictor_vars, synchro_vars)],
  scale.unit = TRUE,
  ncp = 5,
  graph = FALSE,
  quanti.sup = (length(predictor_vars) + 1):(length(predictor_vars) + length(synchro_vars))
)



# ---- Visualization ----

fviz_pca_var(res.pca, repel = TRUE, col.var = "black") +
  ggtitle("Predictor PCA with responses as supplementary variables")






