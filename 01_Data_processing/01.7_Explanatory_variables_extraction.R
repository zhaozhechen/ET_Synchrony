# Author: Zhaozhe Chen
# Date: 2025.11.4

# This code is to extract explanatory variables from remote sensing dataset
# Data sources should be replaced by data at their native spatial scales

library(raster)
library(sf)
# --------- Global --------
# Combined raster from ECOSTRESS results
combined_raster <- readRDS("D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/ECOSTRESS_partial_results/Combined_raster.rds")
# Updated site info
site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

source("01_Data_processing/AMF_processing_functions.R")
Output_path <- "00_Data/"

# ------ Main ---------
# Get site coordinates
coords <- site_info[,c("longitude","latitude")]

# Extract AI values
AI_values <- raster::extract(combined_raster[["AI"]],coords,
                             method = "bilinear")
# Extract canopy height
CH_values <- raster::extract(combined_raster[["CH"]],coords,
                             method = "bilinear")
# Rooting depth
RD_values <- raster::extract(combined_raster[["RD"]],coords,
                             method = "bilinear")
# Sand fraction
TSand_values <- raster::extract(combined_raster[["T_Sand"]],coords,
                                method = "bilinear")
# Put them together into 1 df
predictor_df <- data.frame(site_ID = site_info$site_id,
                           AI = AI_values,
                           CH = CH_values,
                           RD = RD_values,
                           TSand = TSand_values)
# Output this df
write.csv(predictor_df,paste0(Output_path,"perdictor_df.csv"))
