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
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Output path for figures
Output_path <- "03_PN_construction/Results/"
my_color <- brewer.pal(8,"Set2")

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
# Match Site PC scores with coordinates
site_scores <- site_scores %>%
  left_join(Site_info %>%
              select(c(site_id,latitude,longitude)),
            by = c("site_ID" = "site_id"))
# Color coded by PC1
map_PC1 <- syc_map(df = site_scores,
                   varname = "Dim.1",
                   palette_name = "RdYlBu",
                   legend_title = "PC1",
                   g_title = "",
                   color_limits = c(-3,3))
# Color coded by PC2
map_PC2 <- syc_map(df = site_scores,
                   varname = "Dim.2",
                   palette_name = "RdYlBu",
                   legend_title = "PC2",
                   g_title = "",
                   color_limits = c(-3,3))
g <- plot_grid(map_PC1,map_PC2,labels="auto")
print_g(g,"PC_map_syc_metrics_16_psi_ET_FT",10,2.5)

# k-means clustering ======================
# Extract PCA site scores
pc_scores <- site_scores %>% select(Dim.1, Dim.2, Dim.3, Dim.4, Dim.5)
# Elbow method
g_elbow <- fviz_nbclust(pc_scores, kmeans, method = "wss") +
  labs(title = "Elbow Method for k-means")
print_g(g_elbow,"Elbow_syc_metrics_16_psi_ET_FT",4,4)

# Test k numbers from 4-8
set.seed(111) 
for(k in 4:8){
  # k-means
  km_res <- kmeans(pc_scores, centers = k, nstart = 25)
  # Add cluster info to site_scores
  site_scores$cluster <- factor(km_res$cluster)
  
  # Make cluster plot
  g_cluster <- cluster_pt(site_scores)
  # Make a map color coded by cluster
  map_cluster <- syc_map_disc(df = site_scores,
                              varname = "cluster",
                              legend_title = "Cluster",
                              g_title = "")
  # Combine these two plots
  g <- plot_grid(g_cluster,map_cluster,nrow=1,labels = "auto")
  title <- paste0("Cluster",k,"_syc_metrics_16_psi_ET_FT")
  print_g(g,title,10,4)
}





