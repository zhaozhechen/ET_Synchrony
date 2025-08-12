# Author: Zhaozhe Chen
# Date: 2025.8.12

# This code is to smooth LAI using Savitzky-Golay filter
# And calculate start of season (SOS) and end of season (EOS) for each year, at each sites

# --------- Global ----------
library(signal) # the function sgolayfilt is used to apply Savitzky-Golay smoothing filter
library(dplyr)
library(here)
library(zoo) # the function na.approx is used to replace NA by interpolation

# Input path for LAI data
LAI_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/AMF_LAI/MCD15A3H_Lai_500m_"
# Input path to AMF site info
site_info <- read.csv(here("00_data","ameriflux_site_info_update.csv")) %>% select(-X)
# Import ploting functions
source(here("05_Visualization/Plotting_functions.R"))
# parameters to use in the SG filter
# window size
windowsize <- 13
# degree of polynomial
degree_p <- 4
# Colors for plotting
my_color <- brewer.pal(3,"Set2")


# --------- Main ---------
for(i in 1:nrow(site_info)){
  Site_ID <- site_info$site_id[i]
  # Get LAI data
  LAI_df <- read.csv(paste0(LAI_path,Site_ID,".csv"))
  LAI_df <- LAI_df %>%
    mutate(
      # Remove "F" and make it NA
      Lai = as.numeric(if_else(Lai == "F",NA,Lai)),
      # Fill small gaps, fill <=2 consecutive gaps with NA
      LAI_filled = na.approx(Lai,na.rm=FALSE),
      # Apply Savistzky-Golay smoothing filter
      LAI_smoothed = sgolayfilt(LAI_filled,p = degree_p,n = windowsize),
      Date = as.Date(sub("^A","",Date),"%Y%j")
    )
  # Calculate LAI thresholds (50% of max and min LAI) for each year
  GS_df <- LAI_df %>%
    group_by(year(Date)) %>%
    mutate(LAI50 = 0.5*(max(LAI_smoothed)+min(LAI_smoothed))) %>%
    summarise(SOS = Date[which(LAI_smoothed > LAI50)[1]],
              EOS = Date[rev(which(LAI_smoothed > LAI50))[1]]) %>%
    ungroup()
    

  
  
  
  # Calculate SOS and EOS
  
  #Make plots
  g <- ggplot(data=LAI_df,aes(x = Date))+
    geom_line(aes(y=LAI_filled,color="Raw"),size=1)+
    geom_line(aes(y=LAI_smoothed,color="Smoothed"),size=1)+
    my_theme+
    scale_color_manual(values = my_color[1:2],
                       labels = c("Raw","Smoothed"))+
    labs(x="",y="LAI",color="")+
    theme(legend.position = c(0.85,0.9),
          legend.background = element_blank())
  
  # Match with site_info
}


