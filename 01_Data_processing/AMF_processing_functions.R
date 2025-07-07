# Author: Zhaozhe Chen
# Date: 2025.7.6

# These are functions for pre-processing AMF dataset
# 01.1_AMF_processing.R is the main function which outputs a df of required hourly variables
# 01.2_AMF_secondary_variables.R is the main function which XXXX

library(dplyr)
library(ncdf4)

# Get required variable after QC
# For fine resolution (Hourly or Half-hourly), QC = 0 and QC = 1 were kept
# For coarser resolutions (daily through yearly), QC <=0 were removed
# Input include:
# The exact variable name (which should be the sub-variable name for a variable group)
# df: The data frame for processing, normally the raw AMF dataset
# Coarse: whether the temporal resolution is coarse or fine (sub-daily). This should be a logical value (TRUE or FALSE)
Var_QC <- function(varname,df,Coarse = FALSE){
  # Get the variable
  var_values <- df[[varname]]
  # Get the QC name for this variable
  varname_qc <- paste0(varname,"_QC")
  # Get the corresponding QC
  var_QC <- df[[varname_qc]]
  
  if(!varname_qc %in% names(df)){
    warning(paste("Variable",varname,"Does not have a QC file. Returning unmodified values"))
    return(var_values)
  }
  
  if(Coarse){
    # Remove all QC <= 0
    var_values[var_QC <= 0] <- NA
  }else{
    # Keep QC ==0|1
    var_values[var_QC != 0 & var_QC != 1] <- NA
  }
  message(paste("Complete QC for",varname))
  return(var_values)
}

# Since there could be multiple sub-variables for each variable
# This function extracts all sub-variable for the required variable, QC each of them, and get the average of them
# At hourly scale
# Input include:
# key_varname: Part of the variable name, which could be used to uniquely identify the variable
# df: The data frame, which should be the raw AMF dataset before QC
# Coarse: if True, the resolution is daily, if FALSE, the resolution is sub-daily (default: FALSE)
Mean_QC_Var <- function(key_varname,df,Coarse = FALSE){
  # Identify all variable names that starts with key_varname
  varnames <- grep(paste0("^",key_varname),names(df),value=TRUE)
  # Exclude those end with _QC
  varnames <- varnames[!grepl("_QC$",varnames)]
  
  if(length(varnames) == 0){
    stop(paste("No matching variable in this dataset for",key_varname))
  }
  
  # QC for each of the sub-variables
  qc_vars <- lapply(varnames,function(var) Var_QC(var,df,Coarse = Coarse))
  
  # Get means of the sub-variables
  qc_matrix <- do.call(cbind,qc_vars)
  qc_var_mean <- rowMeans(qc_matrix,na.rm = TRUE)
  return(qc_var_mean)
}

# This function is the wrapped function to extract all required variables and QC them
# Input includes:
# var_ls: the list of variables to process
# df: the dataset to be processed, usually the raw AMF dataset, with original variables and their QC variables
Var_QC_all <- function(var_ls,df){
  # Conduct QC for all required variables, and get the mean of them after QC
  results_ls <- lapply(var_ls, function(var){
    Mean_QC_Var(var,df,Coarse = FALSE)
  })

  # Combine these variables into a data frame
  Var_df <- as.data.frame(results_ls)
  names(Var_df) <- var_ls
  Var_df <- cbind(Site_ID = df$site_id,
                  Time = ymd_hms(df$date_time,tz="UTC"),
                  Var_df)
    
  # Set VPD, SWC and TS to 0 if negative
  Var_df$VPD[Var_df$VPD < 0] <- 0
  Var_df$SWC[Var_df$SWC < 0] <- 0
  Var_df$TS[Var_df$TS < 0] <- 0
  return(Var_df)
}

# This function is to standardize time, when 00:00:00 is removed in Time
# Input the dataframe, which should include a Time column as "YYYY-MM-DD hh:mm:ss"
Standardize_time <- function(df){
  # Convert Time to POSIXct Time
  df <- df %>%
    mutate(Time = if_else(
      nchar(Time) == 10, # If it is YYYY-MM-DD
      paste(Time,"00:00:00"),
      Time
    )) %>%
    mutate(Time = ymd_hms(Time,tz="UTC")) %>%
    # Extract hour of the day
    mutate(Hour = hour(Time),
           Date=as.Date(Time)) 
}

# This function finds the nearest grid cell that is not 0 and gets value from that cell
# The input include:
# site_lat: the latitude of the site
# site_lon: the longitude of the site
# grid_value: the values to be matched to, as a matrix
# grid_lat and grid_lon: the corresponding coordinators for the matrix
# max_r: the search radius for values
get_nearest_value <- function(site_lat,site_lon,grid_value,grid_lat,grid_lon,max_r){
  # Find the index for the nearest lat
  lat_idx <- which.min(abs(site_lat - grid_lat))
  # Find the index for the nearest lon
  lon_idx <- which.min(abs(site_lon - grid_lon))
  # Get the value from that cell
  value <- grid_value[lat_idx,lon_idx]
  # Check if this value is NA
  if(!is.na(value) && value !=0){
    return(value)
  }else{
    # If the center cell is 0, search for nearby cells
    for(r in 1:max_r){
     lat_range <- (lat_idx - r):(lat_idx + r)
     lon_range <- (lon_idx - r): (lon_idx + r)
     
     for(i in lat_range){
       for(j in lon_range){
         val <- grid_value[i,j]
         if(!is.na(val) && val !=0){
           return(val)
         }
       }
     }
    }
    # If all nearby values are 0 or NA, return NA
    return(NA)
  }
}

# This function is to extract lon,lat,and values from nc files
# Input include:
# nc_path: the folder path of the nc file
# varname: the key word in the nc file name
extract_nc <- function(nc_path,varname){
  # Get the file path of the target nc file
  nc_file_path <- list.files(nc_path,full.names = TRUE)[grepl(varname,list.files(nc_path))]
  nc <- nc_open(nc_file_path)
  # Get coordinates
  lon <- ncvar_get(nc,"lon")
  lat <- ncvar_get(nc,"lat")
  # Get target variable from the top layer
  var_value <- ncvar_get(nc,varname)[1,,]
  return(list(lon=lon,lat=lat,layer=var_value))
}

# This function calculates soil water potential psi_soil based on soil texture and soil moisture
# Using Clapp and Hornberger equation: psi_soil = psi_ae*(SM/porosity)^(-b)
# Reference: Lauren et al. 2023
# Input include:
# df: the target hourly df of AMF data
# siteinfo: siteinfo including the soil properties
# Site_ID: ID of this site (note: this should be revised later for full AMF dataset)
Cal_psisoil <- function(df,siteinfo,Site_ID){
  # air-entry water potential Unit: kPa
  psi_ae <- siteinfo$psi.ae[siteinfo$site_id == Site_ID]
  # porosity
  porosity <- siteinfo$porosity[siteinfo$site_id == Site_ID]
  # pore-size distribution lamda
  lamda <- siteinfo$lamda[siteinfo$site_id == Site_ID]
  # b = 1/lamda
  b <- 1/lamda
  # SM from AMF data Unit: m3/m3
  SM <- df$SWC/100
  psi_soil <- psi_ae*((SM/porosity)^(-b))
  return(psi_soil)
}
