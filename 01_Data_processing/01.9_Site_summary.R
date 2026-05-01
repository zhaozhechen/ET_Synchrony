# Author: Zhaozhe Chen 
# Date: 2026.4.30
# This code is to summarize site data

# ---- Global ----------
Site_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Data/02_AMF_cleaned/AMF_Hourly/AMF_sites_hourly_update/"

# Get file name
#file_names <- list.files(Site_path,full.names = TRUE)

# Site info file
site_info <- read.csv("00_Data/ameriflux_site_info_update_GS_LAI_filtered.csv")

source("05_Visualization/Plotting_functions.R")

Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Var_summary_plots/"

# Aridity index map
AI_raster <- rast("00_Data/2022_CONUS_AI.tif")
# ----- Main -------
# Sites to test
site_ls <- site_info$site_id

# Initialize vectors to store results
n_year_all <- c()
n_obs_all <- c()
for(arrayid in 1:length(site_ls)){
  # Get file name for each site
  file_name <- paste0(Site_path,"AMF_hourly_",site_ls[arrayid],".csv")
  
  # Read in the df
  site_df <- read.csv(file_name)
  
  site_df$Time <- as.Date(site_df$Time)
  
  # Start time
  t_start <- site_df$Time[1]
  t_end <- tail(site_df$Time,1)
  
  # Get years
  n_year <- as.numeric(difftime(t_end, t_start, units = "days")) / 365
  
  # Get sample size
  n_obs <- nrow(site_df)
  
  # Store these results
  n_year_all <- c(n_year_all,n_year)
  n_obs_all <- c(n_obs_all,n_obs)
  
  message(arrayid)
}

sum_df <- data.frame(
  site_id = site_ls,
  n_year = n_year_all,
  n_obs = n_obs_all)

write.csv(sum_df,paste0(Site_path,"01_data_summary.csv"))

# Make a study site map of 156 AMF sites, size by # of years ------------
site_info$n_year <- sum_df$n_year
g <- ggplot()+
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf,
           fill = "white", colour = NA)
CONUS_outer_4326 <- sf::st_union(CONUS) |>
  sf::st_as_sf() |>
  sf::st_transform(4326)
conus_v_4326 <- terra::vect(CONUS_outer_4326)

base_crop <- terra::crop(AI_raster, terra::ext(conus_v_4326))
base_conus <- terra::mask(base_crop, conus_v_4326)

g <- g +
  tidyterra::geom_spatraster(data = base_conus, alpha = 0.7)
g <- g +
  scale_fill_distiller(
    name = "Aridity Index",
    palette = "RdYlBu",
    direction = 1,
    limits = c(0,2),
    oob = scales::squish,
    na.value = "white",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  )+
  ggnewscale::new_scale_fill()+
  geom_sf(data = sf::st_transform(CONUS, 4326),
          fill = NA, color = "black", alpha = 0.6) +
  coord_sf(crs = sf::st_crs(4326), expand = FALSE)

g <- g +
  geom_point(
    data = site_info,
    aes(x = longitude, y = latitude,
        size = n_year),
    shape = 21,
    fill = "seagreen3",
    color = "black",
    alpha = 0.6
  )+
  scale_size_continuous(
    name = "Number of years",
    breaks = c(5, 10, 15, 20), 
    range = c(2, 8)
  )+
  map_theme +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )
# Output the map
print_g(g,paste0("156_Site_map"),8,5)

