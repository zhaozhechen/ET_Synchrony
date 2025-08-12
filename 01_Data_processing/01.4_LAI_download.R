# Author: Zhaozhe Chen
# Date: 2025.8.11

# This code is to batch download LAI at AMF sites
# Reference: https://modis.ornl.gov/data/modis_webservice.html#apiExamples

# ------- Global --------
library(httr)
library(here)
library(readr)
library(jsonlite)
library(dplyr)
# Base url
url_base <- "https://modis.ornl.gov/rst/api/v1/"
# Product
product <- "MCD15A3H"
network <- "ameriflux"
band <- "Lai_500m"
# Input path to AMF site info
site_info <- read.csv(here("00_data","ameriflux_site_info.csv"))
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/AMF_LAI/"

# -------- Main -------
# Loop over the sites in site_info
for(i in 1:nrow(site_info)){
  # Get Site ID
  Site_ID <- site_info$site_id[i]

  # Get the start date, in the format of AYYYYDOY
  start_d <- paste0("A",site_info$data_start[i],"001")
  end_d <- paste0("A",site_info$data_end[i]+1,"001")

  # Get the full url
  full_url <- paste0(url_base,product,"/",network,"/",Site_ID,"/subsetFiltered?startDate=",
                     start_d,"&endDate=",end_d,"&band=",band)
  r <- GET(full_url,add_headers(Accept = "text/csv"))  
  subset <- read_csv(content(r,as="text"),col_names = FALSE,show_col_types = FALSE)
  # Only keep Date and the central pixel, which should be pixel #145, Column #150
  subset <- subset[,c(3,150)]
  names(subset) <- c("Date","Lai")
    
  # Output this file
  write.csv(subset,paste0(Output_path,product,"_",band,"_",Site_ID,".csv"))
  print(i)
}

