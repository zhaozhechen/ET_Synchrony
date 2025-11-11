# Author: Zhaozhe Chen
# Date: 2025.11.5

# This code is to plot Process Network Chord diagram
# Reference: https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html#directional-relations

# ----- Global -----------
library(dplyr)
library(circlize) # For chord diagram
library(RColorBrewer)

# Input path for Synchrony metrics for 12 pairs
Syc_metrics_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Syc_metrics_12pairs/"

# Updated site info
Site_info <- read.csv("00_Data/ameriflux_site_info_update_GS.csv")

# All variable pairs to consider
var_ls <- c("ET","psi","VPD","TA")
# All seasons
seasons <- c("FT","GS","NGS")

# All variable combinations
var_comb <- expand.grid(from = var_ls,
                        to = var_ls) %>%
  filter(from != to)

# node order around the circle
var_order <- c("ET","psi","VPD","TA")
# Color for the ribbons 
my_colors <- RColorBrewer::brewer.pal(4,"Set3")
cols <- c(ET=my_colors[1], psi=my_colors[2], VPD=my_colors[4], TA=my_colors[3]) 

# Output path for figures
Output_path <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/Results/Hourly_TE_all_sites_server/Results/Chord_diagrams/"

# ------ Main -------
for(i in 1:163){
  Site_ID <- Site_info$site_id[i]
  season <- "GS"
  # This is the variable name for the target synchrony metric
  target_syc_metric_name <- "daily_p_TE"
  
  # Data processing ========================
  # Initialize a vector to store target syc metrics across source-sink pairs
  syc_metrics_pairs <- c()
  # Get target synchrony metric between all source-sink pairs, for the target season
  for(arrayid in 1:nrow(var_comb)){
    # Source variable name
    source_name <- var_comb$from[arrayid]
    # Sink variable name
    sink_name <- var_comb$to[arrayid]
    
    # Read in synchrony metrics df
    syc_df <- read.csv(paste0(Syc_metrics_path,"Syc_metrics_df_",source_name,"_",sink_name,".csv"))
    
    # Extract only the target variable, for this site
    target_syc_metric <- syc_df[[paste0(season,"_",target_syc_metric_name)]][syc_df$site_ID == Site_ID]
    
    syc_metrics_pairs <- c(syc_metrics_pairs,target_syc_metric)
  }
  
  # Add target syc metrics to var_comb
  var_comb$syc_metric <- syc_metrics_pairs
  
  # Making chord diagram ======================
  # Convert var_comb, so that there are two values for the two directions between each source-sink pair
  df_plot <- var_comb %>%
    rename(value_direct1 = syc_metric) %>%
    left_join(var_comb %>%
                rename(value_direct2 = syc_metric,
                       from = to,
                       to = from),
              by = c("from","to")) %>%
    mutate(from_chr = as.character(from),
           to_chr = as.character(to)) %>%
    filter(from_chr < to_chr) %>%
    dplyr::select(-from_chr,-to_chr) %>%
    # Remove NA
    filter(!is.na(value_direct1) & !is.na(value_direct2))
  
  
  # Ribbon color = whichever direction is stronger for that pair
  ribbon_cols <- ifelse(
    df_plot$value_direct1 >= df_plot$value_direct2,
    cols[as.character(df_plot$from)],   # color by 'from'
    cols[as.character(df_plot$to)]      # color by 'to'
  )
  
  
  png(paste0(Output_path,"Chord_diagram_",season,"_",Site_ID,".png"),
      width = 2200,height=2000,res=200)
  circos.clear()
  
  circos.par(start.degree = 90, 
             gap.after = rep(4, length(unique(c(df_plot$from,df_plot$to)))),
             track.height = 0.05)
  
  chordDiagram(
    x = df_plot[,1:4],
    order = var_order,
    grid.col = cols,
    col = ribbon_cols,
    transparency = 0.5,
    directional = 1,
    diffHeight = 0,
    link.sort = TRUE,
    link.largest.ontop = TRUE,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.001)
  )
  
  circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    xlim   <- get.cell.meta.data("xlim")
    ylim   <- get.cell.meta.data("ylim")
    
    # smaller ticks
    circos.axis(h = "top", labels = FALSE, major.tick.length = convert_y(3, "mm"))
    
    # bigger labels; psi as expression
    lab_map <- list(ET = expression(Delta~ET), VPD = expression(Delta~VPD), TA = expression(Delta~"T"), psi = expression(Delta~psi))
    circos.text(mean(xlim), ylim[1] + 50, labels = lab_map[[sector]],
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 2)
  })
  dev.off()
  message(i)
}




