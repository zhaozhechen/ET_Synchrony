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

# Global variables ============
# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")
# All variable combinations
var_comb <- expand.grid(from = var_ls,to = var_ls) %>%
  filter(from != to)

# List of synchrony metric names. Note: if needed, add additional syc metrics here
target_syc_metric_name_ls <- c("daily_agg_TE","daily_p_TE")
# List of their name labels
target_syc_metric_label_ls <- list(bquote(DailyAggTE~"("~psi~"→"~ET~")"),bquote(TE[max]))

# Continuous Predictors to include. Note: these are the initial predictors to include.
pre_var_cont_ls <- c("AVG_SILT","AVG_CLAY","elevation","porosity","CH_GLAD","AI","RD",
                     "avg_air_temp_degC","avg_precip_mm")

# Labels for all predictors
pre_var_cont_label_ls <- list("Silt (%)","Clay (%)","Elevation","Porosity","Canopy height","Aridity index","Rooting depth",
                              bquote("Long-term"~T[air]),"Long-term P")

# Determine metrics to process ===========
# arrayid determines the source and sink variable 
arrayid <- 1
# season includes "GS","NGS","FT". Determines which season to consider
season <- "GS"
# Determine which syc metric to process
syc_id <- 2
# Log transformation of elevation?
log_ele <- FALSE

# Remove highly correlated variables from the predictors (r > 0.7)
# Remove Silt (highly correlated with clay and porosity)
# Remove Long-term P (highly correlated with AI)
remove_id <- c(1,9)

# ------ Main ----------
# Step 1. Process synchrony metrics and match with predictors ============
# log transformation of elevation
if(log_ele == TRUE){
  predictor_df$elevation <- log(predictor_df$elevation)
  pre_var_cont_label_ls[pre_var_cont_label_ls=="Elevation"] <- "log(Elevation)"
}

# Get target synchrony metric name
target_syc_metric_name <- paste0(season,"_",target_syc_metric_name_ls[syc_id])
target_syc_metric_label <- target_syc_metric_label_ls[[syc_id]]

# Extract target synchrony metric values
source_name <- var_comb$from[arrayid]
sink_name <- var_comb$to[arrayid]
# Read in synchrony metrics df between these two pairs
syc_df <- read.csv(paste0(syc_metrics_path,"Syc_metrics_df_",source_name,"_",sink_name,".csv"))
# Extract only the target synchrony metric
syc_df <- syc_df %>%
  select(site_ID,target_syc_metric_name) %>%
  rename(site_id = site_ID)
# Get required predictors and merge with response (synchrony metric)
syc_df <- syc_df %>%
  left_join(predictor_df %>%
              select(c(site_id,pre_var_cont_ls,Koppen_clim_class,IGBP_veg)),
            by="site_id")

# Step 2. Correlation matrix ==============
# Construct correlation matrix
CM <- Hmisc::rcorr(as.matrix(syc_df[c(target_syc_metric_name,pre_var_cont_ls)]),type = "spearman")
# Labels for the variables
my_label <- c(target_syc_metric_label,pre_var_cont_label_ls)
# Plot correlation matrix
g_CM <- plot_CM(CM,my_label)
print_g(g_CM,paste0("CM_",target_syc_metric_name,"_",source_name,"_",sink_name),10,9)

# Plot bivariate correlation between syc metric and predictors
cor_df <- broom::tidy(CM) %>%
  filter(column2 == target_syc_metric_name) %>%
  mutate(p_color = ifelse(p.value < 0.05, "p < 0.05","p > 0.05"))

label_lookup <- setNames(pre_var_cont_label_ls,pre_var_cont_ls)
g_r <- ggplot(data=cor_df,aes(y=reorder(column1,estimate),x=estimate,fill=p_color))+
  geom_col(color="black") +
  my_theme+
  labs(y = "", x = "Spearman ρ",fill="")+
  theme(legend.position = c(0.2,0.8))+
  scale_fill_manual(values = brewer.pal(5,"Set3")[c(4,5)])+
  scale_y_discrete(labels = function(x) label_lookup[x])
print_g(g_r,paste0("Pearson_r_",target_syc_metric_name,"_",source_name,"_",sink_name),6,5)

# Step 3. SLR of syc metric vs predictors ================
# predictor names, excluding highly correlated variables (visual check from above CM results)
pre_var_ls_SLR <- pre_var_cont_ls[-remove_id]
# Corresponding name labels
pre_var_label_ls_SLR <- pre_var_cont_label_ls[-remove_id]

# Initialize a list to store figures
g_scatter_ls <- list()
for(i in 1:length(pre_var_ls_SLR)){
  # Get predictor name
  pred_name <- pre_var_ls_SLR[i]
  # Its label name
  pred_label <- pre_var_label_ls_SLR[[i]]
  
  # Make scatter plots
  g_scatter <- scatter_vars(syc_df,pred_name,target_syc_metric_name,
                            group_name = "IGBP_veg",xtitle = pred_label,ytitle = target_syc_metric_label,
                            my_color = brewer.pal(11,"Set3"))+
    theme(legend.position = "top")
  g_scatter_ls[[i]] <- g_scatter
}

g_scatter <- plot_grid(plotlist = g_scatter_ls,nrow=2,align = "hv")
print_g(g_scatter,paste0("Scatter_",target_syc_metric_name,"_",source_name,"_",sink_name),12,8)

# Step 4. MLR of syc metric vs predictors ================
# Convert categorical data to factor
syc_df$IGBP_veg <- factor(syc_df$IGBP_veg)
syc_df$Koppen_clim_class <- factor(syc_df$Koppen_clim_class)
# Scale continuous predictors
syc_df_scaled <- syc_df
syc_df_scaled[pre_var_ls_SLR] <- scale(syc_df_scaled[pre_var_ls_SLR])

# Formula
MLR_f <- as.formula(paste0(target_syc_metric_name,"~",paste(c(pre_var_ls_SLR),collapse = "+"),"+IGBP_veg+Koppen_clim_class"))
MLR_fit <- lm(MLR_f,data=syc_df_scaled)

# Step 5. PCA ===========================
# Run PCA: predictors = active, synchrony = supplementary quantitative variables
syc_PCA_df <- syc_df[c(paste0(target_syc_metric_name),pre_var_cont_ls)]
# Only keep complete values
syc_PCA_df <- na.omit(syc_PCA_df)

res.pca <- PCA(
  syc_PCA_df,
  scale.unit = TRUE,
  graph = FALSE,
  # Index for quantitative supplementary variables
  quanti.sup = 1
  # Index for qualitative supplementary variables
  #quali.sup = c(8,9)
)

# Find the two PCs with the highest correlations with target syc metrics
PCa <- which(order(res.pca$quanti.sup$coord)==1)
PCb <- which(order(res.pca$quanti.sup$coord)==2)

# Plotting PC1 + PC2
g_PC1 <- fviz_pca_var(res.pca, repel = TRUE, col.var = "black")
# Plotting the PCs where TE loads are higher
g_PC2 <- fviz_pca_var(res.pca,axes = c(PCa,PCb),repel=TRUE,col.var="black")

g_PC <- plot_grid(g_PC1,g_PC2,nrow=1,align="hv")
print_g(g_PC,paste0("PC_",target_syc_metric_name,"_",source_name,"_",sink_name),8,4)






