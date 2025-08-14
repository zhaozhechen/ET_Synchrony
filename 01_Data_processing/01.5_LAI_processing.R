# Author: Zhaozhe Chen
# Date: 2025.8.14

# This code is to smooth LAI using Savitzky-Golay filter
# And calculate start of season (SOS) and end of season (EOS) for each year, at each sites
# Output LAI plots at each site
# Add annual averaged SOS and EOS to the site_info
# Output updated site_info

# --------- Global ----------
library(signal) # the function sgolayfilt is used to apply Savitzky-Golay smoothing filter
library(dplyr)
library(here)
library(zoo) # the function na.approx is used to replace NA by interpolation
library(lubridate)

# Input path for LAI data
LAI_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/AMF_LAI/MCD15A3H_Lai_500m_"
# Input path to AMF site info
site_info <- read.csv(here("00_data","ameriflux_site_info_update.csv")) %>% select(-X)
# Output path for figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites/LAI_plots/"
# Import plotting functions
source(here("05_Visualization/Plotting_functions.R"))
# parameters to use in the SG filter
# window size
windowsize <- 13
# degree of polynomial
degree_p <- 4
# Colors for plotting
my_color <- brewer.pal(3,"Set2")
# LAI cutoff
LAI_th <- 0.5

# --------- Main ---------

# Initialize a vector to store all SOS and EOS
SOS_ls <- c()
EOS_ls <- c()

for(i in 1:nrow(site_info)){
  
  Site_ID <- site_info$site_id[i]
  # Get LAI data
  LAI_df_raw <- read.csv(paste0(LAI_path,Site_ID,".csv"))
  # Preprocessing
  LAI_df <- LAI_df_raw %>%
    mutate(
      # Remove "F" and make it NA
      Lai = as.numeric(if_else(Lai == "F",NA,Lai))
    ) %>%
    # Remove leading NA
    dplyr::filter(cumsum(!is.na(Lai)) > 0) %>%
    mutate(
      Date = as.Date(sub("^A","",Date),"%Y%j"),
      Year = year(Date),
      DOY = yday(Date)
    ) 
  
  # If no available LAI
  if(sum(LAI_df$Lai,na.rm=TRUE)==0){
    SOS <- NA
    EOS <- NA
    message(Site_ID," has all 0 LAI")
  }else{
    # LAI across years =======================
    # Gap fill and Smooth LAI
    LAI_df <- LAI_df %>%
      mutate(
        # Fill small gaps, fill <=2 consecutive gaps with NA
        LAI_filled = na.approx(Lai,na.rm=FALSE),
        # Apply Savistzky-Golay smoothing filter
        LAI_smoothed = sgolayfilt(LAI_filled,p = degree_p,n = windowsize)
      )
    # Calculate LAI thresholds (50% of max and min LAI) for each year
    # And get SOS and EOS
    GS_df <- LAI_df %>%
      group_by(Year) %>%
      mutate(LAI50 = LAI_th*(max(LAI_smoothed,na.rm=TRUE)+min(LAI_smoothed,na.rm=TRUE))) %>%
      summarise(
        LAI50 = first(LAI50),
        SOS = Date[which(LAI_smoothed > LAI50)[1]],
        EOS = Date[rev(which(LAI_smoothed > LAI50))[1]]) %>%
      ungroup()
    GS_df <- na.omit(GS_df)
    # Make LAI TS plot
    g_LAI_years <- plot_LAI_TS(GS_df,LAI_df,my_color)
    
    # Annual LAI cycle ===================
    # Rescale LAI by each year first
    LAI_rescaled_df <- LAI_df %>%
      group_by(Year) %>%
      mutate(
        year_max = max(LAI_smoothed,na.rm=TRUE),
        year_min = min(LAI_smoothed,na.rm=TRUE),
        year_range = year_max - year_min,
        LAI_rescaled = (LAI_smoothed - year_min)/year_range
      )
     
    # Aggregate and smooth LAI
    LAI_annual_df <- LAI_rescaled_df %>%
      # Aggregate across years
      group_by(DOY) %>%
      summarise(LAI_rescaled_annual_mean = mean(LAI_rescaled,na.rm=TRUE)) %>%
      # Remove leading NA
      dplyr::filter(cumsum(!is.na(LAI_rescaled_annual_mean)) > 0) %>%
      mutate(
        # Fill small gaps, fill <=2 consecutive gaps with NA
        LAI_filled = na.approx(LAI_rescaled_annual_mean,na.rm=FALSE),
        # Apply Savistzky-Golay smoothing filter
        LAI_smoothed = sgolayfilt(LAI_filled,p = degree_p,n = windowsize),
        # Add a dummy year
        Year = "2015",
        Date = as.Date(DOY-1,origin = paste0(Year,"-01-01"))
      )
    # Calculate LAI thresholds (50% of max and min LAI) for each year
    # And get SOS and EOS
    GS_annual_df <- LAI_annual_df %>%
      group_by(Year) %>%
      mutate(LAI50 = LAI_th*(max(LAI_smoothed,na.rm=TRUE)+min(LAI_smoothed,na.rm=TRUE))) %>%
      summarise(
        LAI50 = first(LAI50),
        SOS = Date[which(LAI_smoothed > LAI50)[1]],
        EOS = Date[rev(which(LAI_smoothed > LAI50))[1]]) %>%
      ungroup()
    GS_annual_df <- na.omit(GS_annual_df)
    # Make LAI TS plot
    g_LAI_annual <- plot_LAI_TS(GS_annual_df,LAI_annual_df,my_color)+
      scale_x_date(date_breaks = "3 month",date_labels = "%b")+
      theme(legend.position = "none")+
      labs(y = "Rescaled LAI")
    
    # Put two plots together
    g <- plot_grid(g_LAI_years,g_LAI_annual,
                   nrow = 1,
                   align = "h",
                   rel_widths = c(3,1))
    # Output the figure
    print_g(g,paste0("LAI_",Site_ID),
            12,3)
    
    SOS <- as.character(GS_annual_df$SOS)
    EOS <- as.character(GS_annual_df$EOS)
    message("Complete ",Site_ID," (",i," out of ",nrow(site_info),")")  
  }
  
  # Add average SOS and EOS to the site_info df
  SOS_ls <- c(SOS_ls,SOS)
  EOS_ls <- c(EOS_ls,EOS)
}

site_info$SOS <- SOS_ls
site_info$EOS <- EOS_ls
# Add Length of growing season (LGS)
site_info$LGS <- as.numeric(as.Date(site_info$EOS) - as.Date(site_info$SOS))
# Output this updated site_info
write.csv(site_info,"00_Data/ameriflux_site_info_update_GS.csv")






