# Author: Zhaozhe Chen
# Date: 2025.11.4

# This code is to match Clay and Silt fractions with site info
# Data for soil texture/clay/sand/silt is from SSURGO
# There are 22/163 sites having missing mukey, therefore missing clay/silt/sand fractions,replace these sites with group average

# ------- Global ------------
library(terra)
library(sf)
library(foreign) # Read dbf file
library(dplyr)

# Updated site info
site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# Data path for Area- and Depth-Weighted Averages of Selected SSURGO Variables for CONUS
# https://www.sciencebase.gov/catalog/item/631405c8d34e36012efa31ff#stdorder
SSURGO_text <- read.dbf("D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/03_Misc/Text.dbf")

# Data path for MUKEY (SSURGO map unit, 90m resolution)
MUKEY <- rast("D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/03_Misc/MUKEY90m/MapunitRaster_CONUS_90m1.tif")

# Data path for canopy height from GLAD
CH_raster <- rast("D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/03_Misc/Forest_height_2019_NAM.tif")

# ----- Main -------
# Get site coordinates
sites <- site_info %>%
  dplyr::select(id = site_id,longitude,latitude)

site_pts <- st_as_sf(sites,coords = c("longitude","latitude"),crs=4326)

# Extract and match soil texture features ================================
# Extract the mukey at each point
sites$mukey <- terra::extract(MUKEY,site_pts,method="simple")[,2]

# Match clay/silt/sand fractions with sites
sites_texture <- sites %>%
  # Match mukey with texture
  left_join(SSURGO_text,by="mukey") %>%
  rename(site_id = id) %>%
  # Merge with other info
  left_join(site_info,by="site_id") %>%
  # Remove -9999
  mutate(AVG_CLAY = if_else(AVG_CLAY < 0,NA,AVG_CLAY),
         AVG_SAND = if_else(AVG_SAND < 0,NA,AVG_SAND),
         AVG_SILT = if_else(AVG_SILT < 0,NA,AVG_SILT)) %>%
  # Imputation for missing values in AVG_SAND, CLAY, and SILT
  group_by(soil_texture) %>%
  mutate(AVG_SAND = if_else(is.na(AVG_SAND),mean(AVG_SAND,na.rm=TRUE),AVG_SAND),
         AVG_CLAY = if_else(is.na(AVG_CLAY),mean(AVG_CLAY,na.rm=TRUE),AVG_CLAY),
         AVG_SILT = if_else(is.na(AVG_SILT),mean(AVG_SILT,na.rm=TRUE),AVG_SILT)) %>%
  ungroup()

# Extract and match Canopy Height ===================================


