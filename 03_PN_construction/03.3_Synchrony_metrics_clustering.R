# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2025.9.30

# This code is to conduct PCA of multiple synchrony metrics and conduct clustering

# ------ Global ------ 
library(dplyr)
library(factoextra)
library(FactoMineR)
library(ggrepel)

# Path to synchrony metrics df 
# Note: This is a test df for only full time period, psi -> ET
syc_metrics_df <- read.csv("03_PN_construction/Results/Syc_metrics_16_psi_ET_all_sites.csv")

source("05_Visualization/Plotting_functions.R")
# Updated site info
#Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Output path for figures
Output_path <- "03_PN_construction/Results/"
my_color <- brewer.pal(3,"Set2")

# ------- Main ----------
# data processing =========================
metrics_mat <- syc_metrics_df %>%
  select(-c(X,site_ID))

valid_rows <- apply(metrics_mat,1,function(x) !all(is.na(x)))
metrics_mat <- metrics_mat[valid_rows,]

# Get their Site_IDs
Site_ID_valid <- syc_metrics_df$site_ID[valid_rows]

# Scale metrics to have 0 mean and sd 1
metrics_scaled <- scale(metrics_mat)

# Run PCA ==============================
pca_results <- PCA(metrics_scaled,graph=FALSE)

# Get eigenvalues
eig.val <- as.data.frame(get_eigenvalue(pca_results))
eig.val <- eig.val %>%
  mutate(PCs = c(1:nrow(eig.val))) %>%
  # Only keep the first 10 rows
  slice(1:10)
# Make scree plot
g_scree <- scree_pt(eig.val)

# Get site locations
site_scores <- as.data.frame(pca_results$ind$coord) %>%
  mutate(site_ID = Site_ID_valid)
# Get variable loadings
var_loading <- as.data.frame(pca_results$var$coord) %>%
  mutate(var = rownames(pca_results$var$coord))
# Make biplot
g_biplot <- biplot_pt(site_scores,var_loading,5)

# Combine the scree plot and biplot and output
g <- plot_grid(g_scree,g_biplot,nrow=1,labels = "auto")
print_g(g,"PCA_syc_metrics_16_psi_ET_FT",10,5)

# Make maps of PC loadings -------









