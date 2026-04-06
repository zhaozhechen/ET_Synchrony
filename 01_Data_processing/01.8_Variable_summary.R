# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# This code is to make summary plots of input variables


# ------- Global --------------
source("01_Data_processing/AMF_processing_functions.R")
source("05_Visualization/Plotting_functions.R")
# Input path for hourly AMF data
AMF_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/02_AMF_cleaned/AMF_Hourly/AMF_sites_hourly_update/"
# Input path to AMF site info, which also includes soil info
site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")
# Predictor df, with updated AI_gridded
predictor_df <- read.csv("00_Data/perdictor_df_updated.csv")
# Aridity index map
AI_raster <- rast("00_Data/2022_CONUS_AI.tif")
# Input path for raw TE results
TE_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/TE_df/"

# 2 colors for growing season and non-growing season
my_color <- RColorBrewer::brewer.pal(10,"Paired")
# Map
map_color <- RColorBrewer::brewer.pal(n = 11, name = "YlGnBu")

# Sites to process
site_ls <- c("US-Rls","US-xBR")

Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Var_summary_plots/"

# ---- Main ------
# Make maps of Shannon entropy -----------
# Initialize a df to store all H values
H_df <- c()
for(arrayid in 1:nrow(site_info)){
  Site_ID <- site_info$site_id[arrayid]
  # Read in TE list
  TE_ls_site <- readRDS(paste0(TE_path,"TE_df_ls_GS_",Site_ID,".rds"))
  # Get Shannon entropy (H)
  H_var <- sapply(1:12, function(i) TE_ls_site[[i]]$Hy[1])
  # These values correspond to H of ET,ET,ET,psi,psi,psi,VPD,VPD,VPD,TA,TA,TA
  # Get max values for each variable, corresponding to H of ET, psi, VPD, and TA
  H_var <- c(max(H_var[1:3]),max(H_var[4:6]),max(H_var[7:9]),max(H_var[10:12]))
  H_df <- rbind(H_df,H_var)
}
# Rename H_df
H_df <- data.frame(H_df)
names(H_df) <- c("H_ET","H_psi","H_VPD","H_TA")
rownames(H_df) <- NULL
H_df$site_id <- site_info$site_id
# Combine with predictor df
H_df <- merge(predictor_df,H_df,by="site_id")

# Make maps of Shannon entropy ----------
g_map_H_ET <- syc_map(df = H_df,varname = "H_ET",colors = map_color,legend_title = "Shannon Entropy (bit)",
                      g_title = "ET",color_limits = c(0,3),end_marks = "0",base_raster = AI_raster)
g_map_H_psi <- syc_map(df = H_df,varname = "H_psi",colors = map_color,legend_title = "Shannon Entropy (bit)",
                      g_title = "psi",color_limits = c(0,3),end_marks = "0",base_raster = AI_raster)
g_map_H_VPD <- syc_map(df = H_df,varname = "H_VPD",colors = map_color,legend_title = "Shannon Entropy (bit)",
                      g_title = "VPD",color_limits = c(0,3),end_marks = "0",base_raster = AI_raster)
g_map_H_TA <- syc_map(df = H_df,varname = "H_TA",colors = map_color,legend_title = "Shannon Entropy (bit)",
                      g_title = "TA",color_limits = c(0,3),end_marks = "0",base_raster = AI_raster)
# Combine these four maps
g_map_H <- plot_grid(g_map_H_ET,g_map_H_psi,g_map_H_VPD,g_map_H_TA,nrow=2,align="hv",labels="auto")
print_g(g_map_H,"H_maps",14,10)

# Box plots of Shannon entropy across four variables ---------
H_df_long <- H_df %>%
  dplyr::select(site_id, H_ET, H_psi, H_VPD, H_TA) %>%
  tidyr::pivot_longer(
    cols = c(H_ET, H_psi, H_VPD, H_TA),
    names_to = "Variable",
    values_to = "H"
  ) %>%
  dplyr::mutate(
    Variable = factor(
      Variable,
      levels = c("H_ET", "H_psi", "H_VPD", "H_TA")
    )
  )

g_box_H <- ggplot(data = H_df_long,aes(x=Variable,y=H,fill=Variable,
                                       color = Variable))+
  geom_half_violin(alpha=0.5,color=NA)+
  geom_boxplot(width=0.1,color="black",outlier.color = NA)+
  geom_jitter(aes(x = as.numeric(Variable)+0.2),
              position = position_jitter(width=0.1),
              alpha=0.7,size=1)+
  my_theme+
  labs(y = "Shannon entropy (bit)",x="")+
  scale_fill_manual(values = my_color[c(3,1,7,9)])+
  scale_color_manual(values = my_color[c(3,1,7,9)])+
  scale_x_discrete(labels = c(
    "H_ET"  = expression(ET),
    "H_psi" = expression(psi),
    "H_VPD" = expression(VPD),
    "H_TA"  = expression(T[air])
  ))

print_g(g_box_H,"H_comparison",4,3)

# Make TS plots for input variables --------------
# Initialize a list
TS_ls <- list()
for(arrayid in 1:2){
  Site_ID <- site_ls[arrayid]
  
  # Data processing and Calculation of secondary variables ---------------- 
  # Read in hourly data for this site
  AMF_df <- read.csv(paste0(AMF_path,"AMF_hourly_",Site_ID,".csv"))
  # Standardize the time of the df
  AMF_df <- Standardize_time(AMF_df) %>%
    select(-c(X.1,X))
  # Complete time to avoid mismatch in time series alignment
  AMF_df <- AMF_df %>%
    complete(Time = seq(min(Time),max(Time),by="1 hour"))
  
  # Convert psi_soil to log(psi_soil) to reduce skewness
  AMF_df$log10_psi_soil <- log10(AMF_df$psi_soil)
  # Calculate change in delta_log10_psi_soil
  delta_log10_psi_soil <- delta_TS(AMF_df,"log10_psi_soil")
  # Calculate delta_ET
  delta_ET <- delta_TS(AMF_df,"ET")
  # Calculate delta_VPD
  delta_VPD <- delta_TS(AMF_df,"VPD")
  # Calculate delta_TA
  delta_TA <- delta_TS(AMF_df,"TA")
  
  # Put them in a df
  df <- data.frame(
    Time = AMF_df$Time[2:nrow(AMF_df)],
    log10_psi = AMF_df$log10_psi_soil[-1],
    ET = AMF_df$ET[-1],
    VPD = AMF_df$VPD[-1],
    TA = AMF_df$TA[-1],
    delta_log10_psi_soil,
    delta_VPD,
    delta_ET,
    delta_TA)
  
  # Add growing season (GS), the time between SOS and EOS
  # Get SOS for this site
  site_SOS <- format(as.Date(site_info$SOS[site_info$site_id == Site_ID]),"%m-%d")
  site_EOS <- format(as.Date(site_info$EOS[site_info$site_id == Site_ID]),"%m-%d")
  
  df <- df %>%
    mutate(
      # Get month and day
      md = format(Time,"%m-%d"),
      GS = if_else((site_SOS <= md & md <= site_EOS),"GS","Non-GS")
    )
  
  # Make TS and histogram plots for input variables --------------
  # TS for the annual cycle
  g_TS_ET <- TS_annual_v2("delta_ET",df,expression(Delta~ET~(mmday^{-1})),my_color[3])
  g_TS_psi <- TS_annual_v2("delta_log10_psi_soil",df,expression(Delta~psi~"(log scale,"~J~kg^{-1}*")"),my_color[1])+
    theme(legend.position = "none")
  g_TS_VPD <- TS_annual_v2("delta_VPD",df,expression(Delta~VPD~(hPa)),my_color[7])+
    theme(legend.position = "none")
  g_TS_TA <- TS_annual_v2("delta_TA",df,expression(Delta~T[air]~(degree*C)),my_color[9])+
    theme(legend.position = "none")
  TS_ls[[length(TS_ls)+1]] <- g_TS_ET
  TS_ls[[length(TS_ls)+1]] <- g_TS_psi
  TS_ls[[length(TS_ls)+1]] <- g_TS_VPD
  TS_ls[[length(TS_ls)+1]] <- g_TS_TA
}

# Combine all TS plots
g_TS_all <- plot_grid(plotlist = TS_ls[c(1,5,2,6,3,7,4,8)],nrow=4,labels = "auto",align="hv")
print_g(g_TS_all,"TS_var",8,10)


