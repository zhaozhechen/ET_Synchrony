# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# This code is to make plots

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(gghalves)
library(stringr)
library(sf)
library(ggrepel)
library(wesanderson)
library(tidyterra)
library(ggnewscale)

# Make CONUS boundary
# Whole US map
CONUS <- st_read("00_Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]

# Theme for all plots
my_theme <- theme(
  #axis.line=element_line(color="black"),
  panel.background = element_blank(),
  panel.border = element_rect(colour="black",fill=NA),
  legend.key = element_blank(),
  #legend.key.size = unit(6,"cm"),
  #aspect.ratio = 1/1,
  #legend.key.size = unit(0.3,'cm'),
  legend.text = element_text(size=14),
  plot.title = element_text(size=14),
  axis.text = element_text(size=14),
  axis.title = element_text(size=14),
  legend.position = "none",
  strip.background = element_rect(color="black",fill="grey90"),
  strip.text = element_text(face = "bold",size=14),
  plot.margin = margin(15,15,15,15)
)

# Theme for maps
map_theme <- theme(
  #axis.line=element_line(color="black"),
  panel.background = element_blank(),
  panel.border = element_blank(),
  legend.key = element_blank(),
  #legend.key.size = unit(6,"cm"),
  #aspect.ratio = 1/1,
  #legend.key.size = unit(0.3,'cm'),
  legend.text = element_text(size=14),
  plot.title = element_text(size=14),
  legend.title = element_text(size=14),
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "right"
)


# This function is to print pdf and png figure
# Input is the figure g,title,width, and height
print_g <- function(g,title,w,h){
  pdf(paste0(Output_path,"/",title,".pdf"),
      width=w,height=h)
  print(g)
  dev.off()
  png(paste0(Output_path,"/",title,".png"),
      width=w,height=h,units = "in",
      res=600)
  print(g)
  dev.off()
}

# This function is to make full time series TS plot
# Input includes:
# varname: the variable name in the df
# df: the data frame
# y_title: title of y axis
# my_color: color for GS and non-GS, should be a vector of 2
TS_all <- function(varname,df,y_title,my_color){
  g <- ggplot(data=df)+
    #geom_line(color=my_color)+
    geom_segment(aes(x=Time, xend = Time,y = 0,yend = .data[[varname]],
                     color = GS))+
    my_theme+
    scale_color_manual(values = c("GS" = my_color[1],
                                  "Non-GS" = my_color[2]))+
    labs(x="",y=y_title)
  return(g)
}

# This function is to get annual/diurnal cycle of the target variable
# Input includes:
# varname: the variable name in the df
# df: the data frame
# cycle: "Annual" or "Diurnal"
var_cycle <- function(varname,df,cycle){
  # For annual cycle
  if(cycle == "Annual"){
    df_tmp <- df %>%
      mutate(DOY = yday(Time)) %>%
      # Calculate daily mean across the years
      group_by(DOY) %>%
      summarise(
        Time = as.Date(format(first(Time),"2020-%m-%d")),
        mean = mean(.data[[varname]],na.rm=TRUE),
        sd = sd(.data[[varname]],na.rm=TRUE))
  }else if(cycle == "Diurnal"){
    # For diurnal cycle
    df_tmp <- df %>%
      mutate(Hour = hour(Time)) %>%
      # Calculate hourly mean across the day
      group_by(Hour) %>%
      summarise(
        Time = first(Hour),
        mean = mean(.data[[varname]],na.rm=TRUE),
        sd = sd(.data[[varname]],na.rm=TRUE))
  }
  return(df_tmp)
}

# This function is to make annual or diurnal time series TS plots
# Input includes:
# varname: the variable name in the df
# df_cycle: summarized df of annual or diurnal cycle
# cycle: "Annual" or "Diurnal"
# var_to_plot: "mean" or "sd"
TS_cycle <- function(df_cycle,cycle,y_title,var_to_plot){
  df_cycle$Type <- factor(df_cycle$Type,levels=c("Original","Diurnal mean","Diurnal anomaly"))
  g <- ggplot(df_cycle,aes(x=Time,color=Type))+
    geom_line(aes(y=.data[[var_to_plot]]),linewidth=1,alpha=0.7)+
    my_theme+
    scale_color_manual(values = my_color)+
    labs(x = "",y=y_title,color="")
  if(cycle == "Annual"){
    g <- g+
      scale_x_date(date_breaks = "2 month",date_labels = "%b")
  }
  if(var_to_plot == "sd"){
    y_title <- bquote(.(y_title)~sd)
    g <- g+
      labs(y = y_title)
  }
  return(g)
}

# This function plots the annual cycle of the target variable, color coded by GS
# Input includes:
# varname: the variable name in the df
# df: original df
# my_color: color for GS and non-GS, should be a vector of 2
TS_annual <- function(varname,df,y_title,my_color){
  df_tmp <- df %>%
    mutate(DOY = yday(Time)) %>%
    # Calculate daily mean across the years
    group_by(DOY) %>%
    summarise(
      Time = as.Date(format(first(Time),"2020-%m-%d")),
      mean = mean(.data[[varname]],na.rm=TRUE),
      sd = sd(.data[[varname]],na.rm=TRUE),
      GS = first(GS)) %>%
    # Break the line for plotting
    arrange(Time) %>%
    mutate(group_id = rleid(GS))
  g <- ggplot(df_tmp,aes(x=Time,y=mean,color=GS,group=group_id))+
    geom_ribbon(aes(ymin = mean - sd,ymax = mean + sd,fill=GS),color=NA,alpha=0.3)+
    geom_line(size=1)+
    scale_color_manual(values = c("GS" = my_color[1],
                                  "Non-GS" = my_color[2]))+
    scale_fill_manual(values = c("GS" = my_color[1],
                                  "Non-GS" = my_color[2]))+
    my_theme+
    labs(x="",y=y_title,color="",fill="")+
    scale_x_date(date_breaks = "2 month",date_labels = "%b")+
    theme(legend.position = c(0.5,0.9),
          legend.background = element_blank(),
          legend.direction = "horizontal")
  return(g)
}

# This function plots the diurnal cycle of the target variable, color coded by GS
# Input includes:
# varname: the variable name in the df
# df: original df
# my_color: color for GS and non-GS, should be a vector of 2
TS_diurnal <- function(varname,df,y_title,my_color){
  # For diurnal cycle
  df_tmp <- df %>%
    mutate(Hour = hour(Time)) %>%
    # Calculate hourly mean across the day
    group_by(Hour,GS) %>%
    summarise(
      Time = first(Hour),
      mean = mean(.data[[varname]],na.rm=TRUE),
      sd = sd(.data[[varname]],na.rm=TRUE))
  g <- ggplot(df_tmp,aes(x=Time,y=mean,color=GS,fill=GS))+
    geom_ribbon(aes(ymin = mean - sd,ymax = mean + sd),color=NA,alpha=0.3)+
    geom_line(size=1)+
    scale_color_manual(values = c("GS" = my_color[1],
                                  "Non-GS" = my_color[2]))+
    scale_fill_manual(values = c("GS" = my_color[1],
                                 "Non-GS" = my_color[2]))+
    my_theme+
    labs(x="Hour of the day",y=y_title)
  return(g)
}

# This function is to get the distribution of the input variable
# Following the same processing method as for TE implementation, so the distribution is the same as that for TE input
# Input includes:
# var: the values of the variable
# ZFlag: whether zero adjustment is needed
# nbins: # of bins for discretization
# my_color: color for this variable
Hist_var <- function(var,ZFlag,nbins,my_color,lower_qt){
  upper_qt <- 1-lower_qt
  # Find bounds for the variable
  var_bd <- find_bounds(var[var!=0],lower_qt,upper_qt)
  lower_bd <- var_bd[1]
  upper_bd <- var_bd[2]
  var <- na.omit(var)
  # If zero adjustment is needed
  if(ZFlag == TRUE){
    nonzero_idx <- which(var !=0 & !is.na(var))
    # Get non-zero values
    nonzero_values <- var[nonzero_idx]
    # Get histogram info for nonzero values, accounting for outliers
    h <- histogram(nonzero_values,nbins = nbins - 1,lower_bd,upper_bd)
    # Get the breaks
    bin_breaks <- h$breaks
    # Add 0 to the first one
    bin_breaks <- c(0,bin_breaks)
    bin_idx <- ZeroAdjustment(var,nbins,ths=10e-4,lower_bd,upper_bd)
    bin_counts <- table(bin_idx)
  }else{
    h <- histogram(var,nbins,lower_bd,upper_bd)
    bin_breaks <- h$breaks
    bin_counts <- h$counts
  }
  # Get a df for plots
  hist_df <- data.frame(
    bin_left = head(bin_breaks,-1),
    bin_right = tail(bin_breaks,-1),
    count = as.numeric(bin_counts)
  ) %>%
    mutate(bin_center = (bin_left + bin_right)/2,
           bin_width = (bin_right - bin_left))
  # Make a histogram
  g <- ggplot(hist_df,aes(x=bin_center,y=count))+
    geom_col(width=hist_df$bin_width,fill=my_color,color="black")+
    labs(x="",y="Count")+
    ggtitle(paste0("Quantile = ",lower_qt*100,"%"))+
    my_theme
  
  return(g)
}

# This function is to get the histogram bins for the target variable for the desired GS status (GS or Non-GS)
# Following the same processing method as for TE implementation, so the distribution is the same as that for TE input
# Input includes:
# varname: the variable name
# df: the original df
# ZFlag: whether zero adjustment is needed
# nbins: # of bins for discretization
# gs: "GS" or "Non-GS"
Hist_var_GS <- function(varname,df,ZFlag,nbins,gs){
  # Get values for the required GS type
  var <- df[[varname]][df$GS == gs]
  # Find bounds for the variable
  var_bd <- find_bounds(var[var!=0],lower_qt,upper_qt)
  lower_bd <- var_bd[1]
  upper_bd <- var_bd[2]
  var <- na.omit(var)
  # If zero adjustment is needed
  if(ZFlag == TRUE){
    nonzero_idx <- which(var !=0 & !is.na(var))
    # Get non-zero values
    nonzero_values <- var[nonzero_idx]
    # Get histogram info for nonzero values, accounting for outliers
    h <- histogram(nonzero_values,nbins = nbins - 1,lower_bd,upper_bd)
    # Get the breaks
    bin_breaks <- h$breaks
    # Add 0 to the first one
    bin_breaks <- c(0,bin_breaks)
    bin_idx <- ZeroAdjustment(var,nbins,ths=10e-4,lower_bd,upper_bd)
    bin_counts <- table(bin_idx)
  }else{
    h <- histogram(var,nbins,lower_bd,upper_bd)
    bin_breaks <- h$breaks
    bin_counts <- h$counts
  }
  # Get a df for plots
  hist_df <- data.frame(
    bin_left = head(bin_breaks,-1),
    bin_right = tail(bin_breaks,-1),
    count = as.numeric(bin_counts)
  ) %>%
    mutate(bin_center = (bin_left + bin_right)/2,
           bin_width = (bin_right - bin_left))
  return(hist_df)
}

# This function makes Histogram plots for the target variable
# Following the same processing method as for TE implementation, so the distribution is the same as that for TE input
# Input includes:
# varname: the variable name
# df: the original df
# ZFlag: whether zero adjustment is needed
# nbins: # of bins for discretization
# my_color: color for GS and non-GS, should be a vector of 2
Hist_GS_plot <- function(varname,df,x_title,ZFlag,nbins,my_color){
  # Get hist_df for GS and Non-GS
  hist_df_GS <- Hist_var_GS(varname,df,ZFlag,nbins,"GS")
  
  # There could be cases when there is no NGS
  # Skip if it errors:
  hist_df_NGS <- tryCatch({
    Hist_var_GS(varname,df,ZFlag,nbins,"Non-GS")    
  },error=function(e){
    message("No Non-GS for ",Site_ID)
    NULL
  })
  
  # Base plot of GS
  g <- ggplot()+
    geom_col(data=hist_df_GS,
             aes(x=bin_center,y=count),
             width=hist_df_GS$bin_width,color="black",fill=my_color[1],alpha=0.4)+
    my_theme+
    labs(x = x_title,y="count")
  
  # Add NGS only if successful
  if(!is.null(hist_df_NGS)){
    g <- g+
      geom_col(data=hist_df_NGS,
               aes(x=bin_center,y=count),
               width=hist_df_NGS$bin_width,color="black",fill=my_color[2],alpha=0.4)  
  }
  return(g)
}

# This function plots all TS and histogram for the target variable
# Including the full TS, annual cycle, diurnal cycle, and distribution, color coded by how GS is defined
# varname: target variable name (the original name)
# y_title: y title
# df: The target df
# my_color: a vector of three
# ZFlag: whether zero-adjustment should be applied to this variable
# nbins: # of bins for discretization
var_plot_TS_Hist <- function(varname,y_title,df,my_color,ZFlag,nbins){
  # Plot the full time series of the variable, color code by GS
  g_full_TS <- TS_all(varname,df,y_title,my_color)
  # Plot the annual cycle of the variable, color code by GS
  g_annual_TS <- TS_annual(varname,df,y_title,my_color)
  # Plot the diurnal cycle of the variable, color coded by GS
  g_diurnal_TS <- TS_diurnal(varname,df,y_title,my_color)
  # Get the distribution of data
  # If all data is NA, leave it blank
  if(sum(!is.na(df[[varname]])) == 0){
    g_Hist <- ggplot() + 
      xlim(0,5)+
      ylim(0,5)
  }else{
    g_Hist <- Hist_GS_plot(varname,df,x_title = y_title,ZFlag,nbins,my_color)  
  }
  
  g_ls <- list(g_full_TS,g_annual_TS,g_diurnal_TS,g_Hist)
  
  # Combine these four plots
  #              align = "h",nrow=1,
   #              rel_widths = c(1.5,1,1,1))
  return(g_ls)
}

# This function makes all TS and histogram plots for the target variable
# Including the plots of original var, moving diurnal mean, and moving diurnal anomaly
# Also plots the annual cycle and diurnal cycle of the target variable
# Input includes:
# varname: target variable name (the original name)
# y_title: the original y title (no _anomaly or _mean)
# df: The target df
# mycolor: a vector of three
# ZFlag: whether zero-adjustment should be applied to this variable
# nbins: # of bins for discretization
var_plots_all <- function(varname,y_title,df,my_color,ZFlag,nbins){
  # varname of moving window mean for this variable
  varname_mean <- paste0(varname,"_mean")
  # varname of moving window diurnal anomaly
  varname_anomaly <- paste0(varname,"_anomaly")
  # Get titles
  y_title_mean <- substitute(y_title~diurnal~mean,list(y_title=y_title))
  y_title_anomaly <- substitute(y_title~diurnal~anomaly,list(y_title=y_title))
  
  # Plot the full time series of the three variables
  g_original <- TS_all(varname,df,y_title,my_color[1])
  g_mean <- TS_all(varname_mean,df,y_title_mean,my_color[2])
  g_anomaly <- TS_all(varname_anomaly,df,y_title_anomaly,my_color[3])
  
  # Make histogram using all data from each full TS
  # Following the same processing method as for TE implementation, so the distribution is the same as that for TE input
  lower_qt_ls <- c(0.001,0.005,0.01,0.05,0.1)
  # Initialize lists for storing distribution for different quantile conditions
  g_hist_original_ls <- list()
  g_hist_mean_ls <- list()
  g_hist_anomaly_ls <- list()
  
  # Test multiple quantiles for dealing with outliers
  for(i in 1:length(lower_qt_ls)){
    lower_qt <- lower_qt_ls[i]
    upper_qt <- 1-lower_qt
    g_hist_original <- Hist_var(df[[varname]],ZFlag,nbins,my_color[1],lower_qt)
    g_hist_original_ls[[i]] <- g_hist_original
    
    g_hist_mean <- Hist_var(df[[varname_mean]],ZFlag,nbins,my_color[2],lower_qt)
    g_hist_mean_ls[[i]] <- g_hist_mean
    
    g_hist_anomaly <- Hist_var(df[[varname_anomaly]],ZFlag,nbins,my_color[3],lower_qt)
    g_hist_anomaly_ls[[i]] <- g_hist_anomaly
  }
  
  # Annual cycle
  # Make a df for annual cycle data
  df_annual <- rbind(var_cycle(varname,df,"Annual"),
                     var_cycle(varname_mean,df,"Annual"),
                     var_cycle(varname_anomaly,df,"Annual"))
  df_annual$Type <- rep(c("Original",
                          "Diurnal mean",
                          "Diurnal anomaly"),
                        each = nrow(df_annual)/3)
  g_annual <- TS_cycle(df_annual,"Annual",y_title,"mean")
  
  # Diurnal cycle
  # Make a df for diurnal cycle data
  df_diurnal <- rbind(var_cycle(varname,df,"Diurnal"),
                      var_cycle(varname_mean,df,"Diurnal"),
                      var_cycle(varname_anomaly,df,"Diurnal"))
  df_diurnal$Type <- rep(c("Original",
                           "Diurnal mean",
                           "Diurnal anomaly"),
                         each = nrow(df_diurnal)/3)
  g_diurnal_mean <- TS_cycle(df_diurnal,"Diurnal",y_title,"mean")
  g_diurnal_sd <- TS_cycle(df_diurnal,"Diurnal",y_title,"sd")
  
  # Put all figures together
  g1 <- plot_grid(plotlist = g_hist_original_ls,nrow=1)
  g2 <- plot_grid(plotlist = g_hist_mean_ls,nrow=1)
  g3 <- plot_grid(plotlist = g_hist_anomaly_ls,nrow=1)
  g4 <- plot_grid(g_annual,g_diurnal_mean,g_diurnal_sd,nrow=1)
  g_all <- plot_grid(g_original,g1,
                  g_mean,g2,
                  g_anomaly,g3,
                  ncol=1,align="v",axis = "lr")
  g_all <- plot_grid(g_all,g4,nrow=2,rel_heights = c(6,1),align="hv")
  return(g_all)
}

# This function is to plot information metrics vs lag
# Input includes:
# TE_df: the result df from TE
# Type: "TE", "MI","Corr","TEnorm"
# my_color: A three vector color for the plots
TE_lag_plot <- function(TE_df,Type,my_color){
  if(Type == "TE"){
    varname <- "TE"
    y_title <- "TE (bits)"
  }else if(Type == "MI"){
    varname <- "MI"
    y_title <- "MI (bits)"
  }else if(Type == "Corr"){
    varname <- "Corr"
    y_title <- "Correlation"
  }else if(Type == "TEnorm"){
    # Normalize TE by Shannon entropy of the sink
    TE_df <- TE_df %>%
      mutate(TEnorm = TE/Hy * 100,
             TEnormcrit = TEcrit/Hy * 100)
    varname <- "TEnorm"
    y_title <- "Uncertainty reduction (%)"
  }
  # Critical value varname
  varname_crit <- paste0(varname,"crit")
  # Plot of the information metric vs lag
  g <- ggplot(TE_df,aes(x=Lag,y=.data[[varname]]))+
    geom_line(linewidth = 0.8,color=my_color[1])+
    geom_line(aes(y=.data[[varname_crit]]),
              linewidth = 0.8,
              linetype = "dashed",
              color=my_color[3])+
    my_theme+
    labs(x = "Lag (hours)",y=y_title,color="")
  # Get the peak and corresponding lag
  if(Type == "TE"|Type == "MI"|Type == "TEnorm"){
    # This should be the maximum value in significant values
    # Only consider the first 24 hours
    TE_df_1day <- TE_df[1:24,]
    
    TE_df_1day[[varname]][TE_df_1day[[varname]] < TE_df_1day[[varname_crit]]] <- NA
    best_lag <- TE_df_1day$Lag[which.max(TE_df_1day[[varname]])]
    # if not all values not significant
    if(length(best_lag) != 0){
      # Add this best lag to the plot, and annotate it
      g <- g +
        geom_vline(xintercept = best_lag,
                   color=my_color[2],
                   linewidth = 0.8)+
        annotate("text",
                 x=best_lag,
                 y=max(TE_df[[varname]]*1.05,na.rm=TRUE),
                 label = paste0("Best Lag = ",best_lag),
                 size=5,
                 hjust = -0.2,
                 color=my_color[2])      
    }
  }
  return(g)
}

# This function makes plots of information metrics vs lag for different metrics
# Input includes:
# TE_df: the result df from TE
# title: title for the plots
lag_plots_all <- function(TE_df,my_title){
  g_TE <- TE_lag_plot(TE_df,"TE",my_color)+
    ggtitle(my_title)
  g_TEnorm <- TE_lag_plot(TE_df,"TEnorm",my_color)
  g_MI <- TE_lag_plot(TE_df,"MI",my_color)
  g_Corr <- TE_lag_plot(TE_df,"Corr",my_color)
  # Put these four plots together
  g_info <- plot_grid(g_TE,g_TEnorm,g_MI,g_Corr,
                      nrow = 1,
                      align = "hv")
  return(g_info)
}

# This function gets the peak normalized TE value and the corresponding lag
# Note: should consider critical values. Need to update
peak_lag <- function(TE_df){
  # Normalize TE by the Shannon entropy of the sink
  TE_df <- TE_df %>%
    mutate(TEnorm = TE_df$TE/Hy * 100)
  peak_TE <- max(TE_df$TEnorm)
  lag <- TE_df$Lag[which.max(TE_df$TEnorm)]
  return(c(peak_TE,lag))
}

# This function makes LAI plot
# Input includes:
# GS_df: the df summarizes LAI50, SOS, ESO across years
# LAI_df: the df including filled and smoothed LAI
# my_color
plot_LAI_TS <- function(GS_df,LAI_df,my_color){
  # Make blocks df for plotting GS
  gs_shading_df <- GS_df %>%
    transmute(
      xmin = SOS,
      xmax = EOS,
      ymin = -Inf,
      ymax = Inf
    )
  
  # Make df for plotting LAI50
  #LAI50_df <- GS_df %>%
  #  mutate(xmin = as.Date(paste0(Year,"-01-01")),
  #         xmax = as.Date(paste0(Year,"-12-31")),
  #         y=LAI50)
  
  #Make plots
  g <- ggplot(data=LAI_df,aes(x = Date))+
    geom_rect(data=gs_shading_df,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
              fill = my_color[1],inherit.aes = FALSE,alpha=0.3)+
    geom_line(aes(y=LAI_filled,color="Raw"),size=1)+
    geom_line(aes(y=LAI_smoothed,color="Smoothed"),size=1)+
    #geom_segment(data=LAI50_df,aes(x=xmin,xend=xmax,y=y,yend=y,linetype = "LAI50"),
    #             inherit.aes = FALSE)+
    #scale_linetype_manual(name = "",values = c("LAI50" = "dashed"))+
    my_theme+
    scale_color_manual(values = my_color[c(3,2)],
                       labels = c("Raw","Smoothed"))+
    labs(x="",y="LAI",color="")+
    theme(legend.direction = "horizontal",
          legend.position = "top")
    #theme(legend.position = c(0.85,0.8),
    #      legend.title = element_blank(),
    #      legend.background = element_rect(color="black"))
  return(g)
}

# This function is to make Lomb-Scargle Periodogram for the TE vs lag plots
# Input is the LS results, and the color for the peak line
plot_LSP <- function(LSP_results,p_color){
  df <- data.frame(freq = LSP_results$scanned,
                   power = LSP_results$power)
  peak_freq <- LSP_results$peak.at[1]
  peak_period <- round(1/peak_freq,2)
  peak_power <- LSP_results$peak
  
  g <- ggplot(df,aes(x=freq,y=power))+
    geom_line(linewidth = 0.8)+
    geom_hline(yintercept = LSP_results$sig.level,linetype = "dashed")+
    geom_vline(xintercept = peak_freq,color=p_color,linewidth=0.8)+
    labs(x = "Frequency",y="Normalized Power")+
    my_theme+
    scale_x_continuous(n.breaks=3)
  
  if(peak_power > LSP_results$sig.level){
    g <- g +
      geom_text(aes(x=peak_freq,y=peak_power - 0.02),
                label = paste0("P = ",peak_period,"h"),
                color=p_color,size=5,hjust=-0.2)
  }
  
  return(g)
}

# This function is to make normalized TE plot vs lag, with synchrony metrics annotated
# Input include:
# TE_df_tmp: TE_df with normalized TE and normalized TE critical value
# p_lag: peak lag in hours
# p_TE: peak TE value
# mem: memory in hours
TE_norm_lag_plot <- function(TE_df_tmp,p_lag,p_TE,mem,m_color){
  g_TE <- ggplot(TE_df_tmp,aes(x=Lag,y=TE_norm))+
    geom_line(linewidth=0.8)+
    labs(x="Lag (h)",y="Uncertainty reduction (%)")+
    geom_line(aes(y=TEcrit_norm),
              linewidth = 0.8,
              linetype = "dashed",
              color=m_color[2])+
    # Add peak lag
    geom_vline(xintercept = p_lag,
               linewidth = 0.8,
               color=m_color[5])+
    # Annotate peak lag
    annotate("text",
             x=p_lag,
             y=max(TE_df_tmp$TE_norm*1.05,na.rm=TRUE),
             label = paste0("Best Lag = ",p_lag,"h"),
             size=5,
             hjust = -0.2,
             color=my_color[5])+
    # Annotate peak TE value
    annotate("text",
             x=p_lag,
             y=max(TE_df_tmp$TE_norm*0.95,na.rm=TRUE),
             label = paste0("Peak daily TE = ",round(p_TE,2),"%"),
             size=5,
             hjust = -0.1,
             color=my_color[1])+
    # Add memory
    annotate("segment",
             x=p_lag,xend=p_lag + mem,
             y= mean(TE_df_tmp$TEcrit_norm, na.rm=TRUE)*1.05,
             arrow = arrow(ends = "both",type = "closed",length = unit(10,"pt")),
             color = m_color[3])+
    my_theme
  
  if(mem == max_lag){
    g_TE <- g_TE +
      annotate("text",
               x=p_lag,
               y=mean(TE_df_tmp$TEcrit_norm, na.rm=TRUE)*1.2,
               label = paste0("Memory ≥ ",mem,"h"),
               size=5,
               hjust = -0.1,
               color=my_color[3])
  }else{
    g_TE <- g_TE +
      annotate("text",
               x=p_lag,
               y=mean(TE_df_tmp$TEcrit_norm, na.rm=TRUE)*1.2,
               label = paste0("Memory = ",mem,"h"),
               size=5,
               hjust = -0.1,
               color=my_color[3])
  }
  return(g_TE)
}

# Make maps color coded by the target variable 
# This function if to make map of target synchrony metrics
# df: the input df including lat and long, and the target var
# varname: target variable name
# palette_name: palette name to use
# legend_title: legend title
# g_title: Title
# color_limits: need an input for the global range for three seasons, so maps of different
# seasons share the same color bar
# shape_var determines whether the points should be plotted using shape group
syc_map <- function(df, varname,
                    palette_name = NULL,
                    colors = NULL,
                    legend_title = "",
                    g_title = "",
                    color_limits = NULL,
                    direction = -1,
                    n_breaks = 5,
                    end_marks = c("left", "right"),   # choose: "left", "right", or both, or character(0)
                    left_mark = "\u2264",             # "≤"
                    right_mark = "\u2265",             # "≥"
                    shape_var = NULL,
                    shape_values = c(21,23),
                    base_raster = NULL,
                    base_name = "Aridity index",
                    base_alpha=0.7,
                    base_limits = c(0,2),
                    base_colors = NULL,
                    base_n_breaks = 5,
                    base_palette_name = "RdYlBu",
                    base_direction = 1,
                    size_var = NULL,
                    size_range = c(2,6),
                    size_title = "Total TE"
) {
  
  # ---- checks ----
  if (is.null(palette_name) && is.null(colors)) {
    stop("Provide either `palette_name` or `colors`.")
  }
  if (!is.null(palette_name) && !is.null(colors)) {
    stop("Provide only one of `palette_name` or `colors`, not both.")
  }
  if (is.null(color_limits) || length(color_limits) != 2) {
    stop("`color_limits` must be length 2: c(min, max).")
  }
  end_marks <- intersect(end_marks, c("left", "right"))
  
  # ---- base plot ----
  g <- ggplot()+ annotate(
    "rect",
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    fill = "white", colour = NA
  )
  
  # ---- NEW: local CONUS outer boundary in EPSG:4326 (does NOT change global CONUS) ----
  CONUS_outer_4326 <- sf::st_union(CONUS) |>
    sf::st_as_sf() |>
    sf::st_transform(4326)
  
  conus_v_4326 <- terra::vect(CONUS_outer_4326)
  
  # ---- NEW: add AI raster background (cropped/masked in EPSG:4326) ----
  if (!is.null(base_raster)) {
    
    # base_raster is already EPSG:4326 from your printout
    base_crop <- terra::crop(base_raster, terra::ext(conus_v_4326))
    base_conus <- terra::mask(base_crop, conus_v_4326)
    
    g <- g + tidyterra::geom_spatraster(data = base_conus, alpha = base_alpha)
    
    base_scale_args <- list(
      name   = base_name,
      limits = base_limits,
      oob    = scales::squish
    )
    if (!is.null(base_limits) && length(base_limits) == 2) {
      base_scale_args$breaks <- seq(base_limits[1], base_limits[2], length.out = base_n_breaks)
    }
    
    # --- NEW: RdYlBu (blue = high) for AI ---
    base_scale_args$guide <- guide_colorbar(
      ticks = TRUE,
      frame.colour = "black",
      ticks.colour = "black"
    )
    
    g <- g + do.call(
      scale_fill_distiller,
      c(list(
        palette   = base_palette_name,
        direction = base_direction,
        na.value  = "white"
      ), base_scale_args)
    )
    
    g <- g + ggnewscale::new_scale_fill()

  }
  
  # ---- draw CONUS outlines on top, but DON'T fill (so raster stays visible) ----
  g <- g +
    geom_sf(data = sf::st_transform(CONUS, 4326),
            fill = NA, color = "black", alpha = 0.6) +
    coord_sf(crs = sf::st_crs(4326), expand = FALSE)
  
  # Add points -----
  if(!is.null(shape_var)){
    
    if(!is.null(size_var)){
      g <- g +
        geom_point(
          data = df,
          aes(x = longitude.x, y = latitude.x,
              fill = .data[[varname]],
              shape = .data[[shape_var]],
              size = .data[[size_var]]),
          alpha = 0.8, color = "black"
        ) +
        scale_shape_manual(values = shape_values) +
        scale_size_continuous(
          name = size_title,
          range = size_range
        )
    } else {
      g <- g +
        geom_point(
          data = df,
          aes(x = longitude.x, y = latitude.x,
              fill = .data[[varname]],
              shape = .data[[shape_var]]),
          size = 3.5, alpha = 0.8, color = "black"
        ) +
        scale_shape_manual(values = shape_values)
    }
    
  } else {
    
    if(!is.null(size_var)){
      g <- g +
        geom_point(
          data = df,
          aes(x = longitude.x, y = latitude.x,
              fill = .data[[varname]],
              size = .data[[size_var]]),
          alpha = 0.8, shape = 21, color = "black"
        ) +
        scale_size_continuous(
          name = size_title,
          range = size_range
        )
    } else {
      g <- g +
        geom_point(
          data = df,
          aes(x = longitude.x, y = latitude.x,
              fill = .data[[varname]]),
          size = 3.5, alpha = 0.8, shape = 21, color = "black"
        )
    }
  }
  
  g <- g +
    labs(fill = legend_title) +
    ggtitle(g_title) +
    map_theme
  
  # ---- label function ----
  label_fun <- function(x) {
    labs <- as.character(x)
    if (length(labs) >= 2) {
      if ("left" %in% end_marks)  labs[1] <- paste0(left_mark, " ", labs[1])
      if ("right" %in% end_marks) labs[length(labs)] <- paste0(right_mark, " ", labs[length(labs)])
    }
    labs
  }
  
  # ---- common scale args ----
  scale_args <- list(
    limits = color_limits,
    oob = scales::squish,
    breaks = seq(color_limits[1], color_limits[2], length.out = n_breaks),
    labels = label_fun,
    guide = guide_colorbar(
      ticks = TRUE,
      frame.colour = "black",
      ticks.colour = "black"
    )
  )
  
  # ---- choose scale ----
  if (!is.null(colors)) {
    g <- g + do.call(scale_fill_gradientn, c(list(colours = colors), scale_args))
  } else {
    g <- g + do.call(scale_fill_distiller, c(list(palette = palette_name, direction = direction), scale_args))
  }
  
  
  return(g)
}

# This is a simlar function to make maps, but for discrete values
syc_map_disc <- function(df,
                         regime_var = "regime",
                         regime_colors,
                         legend_title = "",
                         g_title = "",
                         base_raster = NULL,
                         base_name = "Aridity index",
                         base_alpha = 0.7,
                         base_limits = c(0,2),
                         base_palette_name = "RdYlBu",
                         base_direction = 1) {
  
  # ---- base plot ----
  g <- ggplot() +
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = -Inf, ymax = Inf,
             fill = "white", colour = NA)
  
  # ---- CONUS boundary ----
  CONUS_outer_4326 <- sf::st_union(CONUS) |>
    sf::st_as_sf() |>
    sf::st_transform(4326)
  
  conus_v_4326 <- terra::vect(CONUS_outer_4326)
  
  # ---- raster background ----
  if (!is.null(base_raster)) {
    
    base_crop <- terra::crop(base_raster, terra::ext(conus_v_4326))
    base_conus <- terra::mask(base_crop, conus_v_4326)
    
    g <- g +
      tidyterra::geom_spatraster(data = base_conus, alpha = base_alpha)
    
    g <- g +
      scale_fill_distiller(
        name = base_name,
        palette = base_palette_name,
        direction = base_direction,
        limits = base_limits,
        oob = scales::squish,
        na.value = "white",
        guide = guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          title.hjust = 0.5
        )
      )
    
    g <- g + ggnewscale::new_scale_fill()
  }
  
  # ---- CONUS outline ----
  g <- g +
    geom_sf(data = sf::st_transform(CONUS, 4326),
            fill = NA, color = "black", alpha = 0.6) +
    coord_sf(crs = sf::st_crs(4326), expand = FALSE)
  
  # ---- points ----
  g <- g +
    geom_point(
      data = df,
      aes(x = longitude.x, y = latitude.x,
          fill = .data[[regime_var]]),
      size = 3.5,
      shape = 21,
      color = "black",
      alpha = 0.9
    )
  
  # ---- categorical colors ----
  g <- g +
    scale_fill_manual(
      values = regime_colors,
      name = legend_title,
      guide = guide_legend(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5
      )
    )
  
  # ---- theme ----
  g <- g +
    ggtitle(g_title) +
    map_theme +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "horizontal",
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  
  return(g)
}

# This function is to get title and legend from the variable name
get_title<- function(varname){
  # Extract title name
  prefix <- str_extract(varname,"^[^_]+(?:_[^_]+)*?(?=_[^_]+_to_)")
  source_name <- str_extract(varname, "(?<=_)[^_]+(?=_to_)")
  sink_name <- str_extract(varname, "(?<=_to_)[^_]+$")
  if(prefix == "p_TE"){
    prefix_title <- "Peak daily TE"
    legend_title <- "%"
  }else if(prefix == "p_lag"){
    prefix_title <- "Best lag"
    legend_title <- "hour"
  }else if(prefix == "mem"){
    prefix_title <- "Memory"
    legend_title <- "hour"
  }else if(prefix == "period"){
    prefix_title <- "Period"
    legend_title <- "hour"
  }else if(prefix == "agg_TE"){
    prefix_title <- "Aggregated daily TE"
    legend_title <- "%"
  }
  # Get title
  title <- bquote(.(prefix_title)~"("~Delta~.(as.name(source_name))~"\u2192"~Delta~.(as.name(sink_name))~")")
  return(list(legend_title,title))
}


# This function makes maps of target variable for three periods: Full TS, GS, and NGS
# df: Input is the full df, including three periods, lat,long,and target variable
# varname: target variable name
# palette_name: palette name to use
season3_syc_map <- function(df,varname,palette_name,legend_title){
  # Get color range
  color_limits <- range(df[[varname]],na.rm=TRUE)
  
  # Full Time
  FT_map <- syc_map(df[df$GS == "FT",],
                    varname,palette_name,
                    legend_title = legend_title,
                    g_title = "Full Time Series (FT)",
                    color_limits = color_limits)
  # GS
  GS_map <- syc_map(df[df$GS == "GS",],
                    varname,palette_name,
                    legend_title = legend_title,
                    g_title = "Growing Season (GS)",
                    color_limits = color_limits)
  # NGS
  NGS_map <- syc_map(df[df$GS == "NGS",],
                     varname,palette_name,
                     legend_title = legend_title,
                     g_title = "Non Growing Season (NGS)",
                     color_limits = color_limits)
 
  return(list(FT_map,GS_map,NGS_map))
}

# This function is to make histogram of target variable
plot_syc_hist <- function(df,varname,my_color,x_title,g_title){
  g <- ggplot(df,aes(x=.data[[varname]]))+
    geom_histogram(bins=10,
                   color="black",
                   fill=my_color)+
    my_theme+
    labs(x=x_title)+
    ggtitle(g_title)
  return(g)
}

# This function is to compare synchrony metrics across groups
# Input include:
# df: input df containing the target variables
# varname: synchrony metric to be compared
# group_name: groups to be compared
# y_title
# my_color: colors for the group
Hist_Syc <- function(df,varname,group_name,y_title,my_color,g_title){
  df <- df[!is.na(df[[varname]]),]
  # Filter out groups with fewer than 2 data points
  df <- df %>%
    group_by(.data[[group_name]]) %>%
    dplyr::filter(n() >= 2) %>%
    ungroup()
  
  # Reset factor levels to match the filtered data
  df[[group_name]] <- factor(df[[group_name]])
  
  g <- ggplot(data=df,aes(x=.data[[group_name]],
                          y=.data[[varname]]))+
    geom_half_violin(fill=my_color,alpha=0.5,color=NA)+
    geom_boxplot(width=0.1,color="black",outlier.color = NA)+
    geom_jitter(aes(x = as.numeric(.data[[group_name]])+0.2),
                position = position_jitter(width=0.1),
                alpha=0.7,
                color=my_color)+
    my_theme+
    labs(x="",y=y_title)+
    coord_flip()+
    ggtitle(g_title)
  return(g)
}

# This function is to make distribution and group comparison plots 
# for the selected time period for the target metric
# df: full df of three periods
# GS: FT, GS, or NGS
# g_title: title for the first plol, could be "Full Time Series", "Growing season", "Non Growing Season"
# my_color: color for this season
# legend_title: legend for the x axis
Hist_season1 <- function(df,GS,g_title,my_color,legend_title){
  df <- Syc_metrics_df[Syc_metrics_df$GS == GS,]
  
  # Make histogram of target variable
  g_hist <- plot_syc_hist(df,varname,my_color,legend_title,g_title)
  
  # Across AI
  g_AI <- Hist_Syc(df,varname,"AI_level",legend_title,my_color,"AI level")
  
  # Across Soil Type
  g_Soil <- Hist_Syc(df,varname,"Soil_type",legend_title,my_color,"Soil type")
  
  # Across Climate
  g_Clm <- Hist_Syc(df,varname,"Koppen_clim_class",legend_title,my_color,"Climate")
  
  # Across IGBP
  g_IGBP <- Hist_Syc(df,varname,"IGBP_veg",legend_title,my_color,"IGBP")  
  
  return(list(g_hist,g_AI,g_Soil,g_Clm,g_IGBP))
}

# This function is to make summary plots of each target metric
# Input include:
# Syc_metrics_df: the df including all synchrony metrics and required explanatory variables
# varname: target variable to plot
# palette_name: name of the palette for making maps
plot_syc_all <- function(Syc_metrics_df,varname,palette_name){
  # Get legend and title
  legend_title <- get_title(varname)[[1]]
  g_title <- get_title(varname)[[2]]
  
  # Make maps of the target variable for FT, GS, and NGS
  map_seasons <- season3_syc_map(Syc_metrics_df,varname,palette_name,legend_title)
  
  # Compare target metrics across three seasons
  g_syc_seasons <- ggplot(data=Syc_metrics_df,
                          aes(x=GS,y=.data[[varname]],
                              color=GS,fill=GS))+
    geom_half_violin(alpha=0.5,color=NA)+
    geom_boxplot(width=0.1,color="black",outlier.color = NA)+
    geom_jitter(aes(x = as.numeric(GS)+0.2),
                position = position_jitter(width=0.1),
                alpha=0.7)+
    my_theme+
    scale_fill_manual(values = season_color[c(3,1,2)])+
    scale_color_manual(values = season_color[c(3,1,2)])+
    labs(x="",y=legend_title)
  
  # Combine three maps + metrics across seasons
  g_row1 <- plot_grid(map_seasons[[1]],map_seasons[[2]],map_seasons[[3]],g_syc_seasons,nrow=1,
                      labels = "auto")
  
  # Add one overall title
  g_row1 <- plot_grid(
    ggdraw() + draw_label(
      g_title,fontface = "bold",size=16,hjust=0.5
    ),
    g_row1,
    ncol=1,
    rel_heights = c(0.08,1)
  )
  
  # FT histograms
  FT_hist_all <- Hist_season1(Syc_metrics_df,"FT","Full Time Series",season_color[3],legend_title)
  # GS histograms
  GS_hist_all <- Hist_season1(Syc_metrics_df,"GS","Growing Season",season_color[1],legend_title)
  # NGS histograms
  NGS_hist_all <- Hist_season1(Syc_metrics_df,"NGS","Non Growing Season",season_color[2],legend_title)
  # Combine these 15 plots
  g_row2 <- plot_grid(plotlist = c(FT_hist_all,GS_hist_all,NGS_hist_all),
                      nrow=3,
                      labels = letters[5:19])
  
  # Combine row1 and row2
  g_all <- plot_grid(g_row1,g_row2,nrow=2,
                     rel_heights = c(1,3))
  # Output this figure
  print_g(g_all,paste0(varname,"_all"),
          20,16)
}


# This function adds p-value to the Histgram
Hist_Syc_p_value <- function(df,varname,group_name,y_title,my_color){
  # Make plot
  g <- Hist_Syc(df,varname,group_name,y_title,my_color)
  # Get p value
  p_value <- syc_compare(df,varname,group_name)
  # Add p value to the plot
  p_value_text <- paste0("p = ",p_value)
  g <- ggdraw()+
    draw_plot(g,x=0,y=0,width=1,height=1)+
    draw_label(p_value_text,x=0.98,y=0.98,hjust=1,vjust=1)
  return(g)
}

# This function make scatter plots to compare two variables
scatter_vars <- function(df,varname1,varname2,group_name,xtitle,ytitle,my_color){
  # Formula for the line
  f <- as.formula(paste0(varname2,"~",varname1))
  # Make a line
  lm <- lm(f,data=df)
  # Get parameters
  k <- lm$coefficients[2]
  b <- lm$coefficients[1]
  #p <- signif(summary(lm)$coefficients[2,4],2)
  p <- formatC(summary(lm)$coefficients[2,4], format = "e", digits = 2)
  R2 <- round(summary(lm)$r.squared,3)
  
  g <- ggplot(data=df,aes(x=.data[[varname1]],y=.data[[varname2]],color=.data[[group_name]]))+
    geom_point(size=2,alpha=0.7)+
    my_theme+
    labs(x=xtitle,y=ytitle,color="")+
    geom_abline(slope = k,intercept = b)+
    annotate(
      geom = "text",
      x=Inf,y=-Inf,
      label = as.expression(bquote(italic(p) == .(p))),
      hjust=1.2,vjust=-9,
      size=5)+
    annotate(
      geom = "text",
      x=Inf,y=-Inf,
      label = as.expression(bquote(R^2 == .(R2))),
      hjust=1.2,vjust=-10,
      size=5)+
    theme(legend.position = c(0.15,0.85),
          legend.title = element_blank(),
          legend.background = element_blank())+
    scale_color_manual(values=my_color)
  return(g)
}

# Dummy function to wrap a g_matrix into a ggplot object
grob_to_ggplot <- function(g_matrix) {
  grob_obj <- grid::grid.grabExpr(print(g_matrix))
  g <- ggplot() + annotation_custom(grob_obj)
  return(g)
}

# This function makes scree plot from pca results
# Input is the eigen value df
scree_pt <- function(eig.val){
  g_scree <- ggplot(data=eig.val,aes(x=PCs,y=variance.percent))+
    geom_bar(stat="identity",
             fill="#2e7ebb",
             alpha=0.8,
             color="black",
             linewidth=0.5)+
    geom_line(color="black",linewidth=0.8)+
    geom_point(size=2)+
    my_theme+
    labs(x="PCs",y = "Explained Variances (%)")+
    scale_x_continuous(breaks = c(1:10))+
    geom_text(aes(label = paste0(round(variance.percent,1),"%")),
              size=4,
              hjust=0.01,vjust=-0.5)+
    ylim(c(0,45))+
    theme(aspect.ratio = 1/1)
  return(g_scree)
}

# This function makes biplot from pca results
# Input include:
# The site location df
# The variable loading df
# A factor to shrink arrows
biplot_pt <- function(site_scores,var_loading,arrow_factor){
  var_loading <- var_loading %>%
    mutate(PC1 = Dim.1 * arrow_factor,
           PC2 = Dim.2 * arrow_factor)
  
  # Make biplot 
  g_biplot <- ggplot()+
    # Arrows
    geom_segment(data = var_loading,aes(x=0,y=0,xend=PC1,yend=PC2),
                 arrow = arrow(length = unit(0.25,"cm")),
                 color=my_color[2],linewidth = 0.8)+
    # Sites
    geom_point(data = site_scores,aes(x = Dim.1,y=Dim.2),
               color=my_color[1],size=2)+
    # Labels for variables
    geom_text_repel(data = var_loading,
                    aes(x = PC1, y = PC2, label = var),
                    color = my_color[2], size = 5, 
                    box.padding = 0.5,
                    point.padding = 0.5,
                    segment.color = "grey50") +
    geom_hline(yintercept = 0,linetype = "dashed",linewidth = 0.8)+
    geom_vline(xintercept = 0,linetype = "dashed",linewidth = 0.8)+
    my_theme+
    labs(x = paste0("PC1 (", round(eig.val$variance.percent[1],1), "%)"),
         y = paste0("PC2 (", round(eig.val$variance.percent[2],1), "%)")) +
    theme(aspect.ratio = 1)  
  return(g_biplot)
}

# This function makes k-means clustering plot
# Input includes site_score with PC1,PC2, and clusters
cluster_pt <- function(site_scores){
  g_cluster <- ggplot(site_scores, aes(x = Dim.1, y = Dim.2, color = cluster)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = my_color) +
    my_theme +
    labs(x = paste0("PC1 (", round(eig.val$variance.percent[1],1), "%)"),
         y = paste0("PC2 (", round(eig.val$variance.percent[2],1), "%)"),
         color = "Cluster",
         title = "K-means Clusters in PCA Space")+
    theme(aspect.ratio = 1) 
  return(g_cluster)
}

# This function makes scatter plots of syc metrics vs predictors, color coded by season
# Input include merged syc_df, with predictors merged
# the response and predictor variable names
syc_scatter <- function(df,res_varname,pre_varname){
  # Reshape df to long format
  df_long <- df %>%
    # Only keep required variables
    dplyr::select(site_ID,all_of(pre_varname),contains(res_varname)) %>%
    pivot_longer(
      cols = matches(paste0("^(FT|GS|NGS)_",res_varname,"$")),
      names_to = "Season",
      values_to = "Response"
    ) %>%
    mutate(
      # Clean up season names
      Season = sub(paste0("_", res_varname), "", Season),
      Season = factor(Season, levels = c("FT", "GS", "NGS"))
    )
  
  # Calculate R2 for each season
  R2_df <- df_long %>%
    group_by(Season) %>%
    do({
      mod <- lm(Response ~ .data[[pre_varname]], data = .)
      s <- summary(mod)
      tibble(
        R2   = s$r.squared,
        pval = coef(s)[2, 4]  # p-value for slope
      )
    }) %>%
    ungroup() %>%
    mutate(
      # format p-values in ggsignif style
      p_label = case_when(
        pval < 0.001 ~ "***",
        pval < 0.01  ~ "**",
        pval < 0.05  ~ "*",
        TRUE         ~ "ns"
      ),
      label = paste0(Season, ": R² = ", round(R2, 2),
                     " (", p_label, ")")
    )
  
  # Get location for labels
  x_max <- max(df_long[[pre_varname]],na.rm=TRUE)
  y_range <- range(df_long$Response,na.rm=TRUE)
  y_span <- y_range[2] - y_range[1]
  y_top <- y_range[2]
  y_offset <- 0.15 * y_span
  
  # Make scatter plot
  g <- ggplot(data=df_long,aes(x=.data[[pre_varname]],
                               y=Response,color=Season,fill=Season))+
    geom_point(size=4,shape=21,alpha=0.8,color="black")+
    # Add fitted line
    geom_smooth(method = "lm",se = TRUE,linetype = "solid",alpha = 0.25)+
    scale_color_manual(values = season_color)+
    scale_fill_manual(values = season_color)+
    geom_text(
      data = R2_df,
      aes(x = x_max,
          y = y_top + 0.6 * y_span - (as.numeric(factor(Season))) * y_offset,
          label = label, color = Season),
      hjust = 1, vjust = 1, size = 5, show.legend = FALSE
    ) +
    my_theme+
    labs(x = pre_varname,y = res_varname,color="")
  
  return(g)
}

# This function is to make heatmap for adjusted R2 or k
# Input include name of variable to be filled
# season: FT, GS, NGS
# palette name for filling
heatmap_plot <- function(fill_name,season,palette_name,direction){
  g <- ggplot(results_all %>% filter(Season == season),
              aes(x = predictor,y = response,
                  fill= if_else(pval < 0.05,.data[[fill_name]],NA)))+
    geom_tile(color="white",linewidth = 0.3)+
    scale_fill_distiller(
      palette = palette_name,
      na.value = "grey",
      direction = direction,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black"
      )
    )+
    facet_wrap(~ pair,ncol=3)+
    my_theme+
    labs(x = "",y="",fill=fill_name)+
    theme(legend.position = "right")+
    ggtitle(season)  
  return(g)
}

# This function is to make density curves for target variable
# Input include:
# df: input df
# varname: target variable name
# group_varname: varname used for color
# color_pallete
# x_intercept: intercept for plotting vertical line
# x_title
plot_pdf <- function(df,varname,group_varname,my_pallete,x_intercept,x_title){
  g <- ggplot(df,aes(x=.data[[varname]],fill=.data[[group_varname]],color=.data[[group_varname]]))+
    geom_density(alpha = 0.4,linewidth = 1)+
    geom_vline(xintercept = x_intercept,linetype = "dashed",linewidth = 1)+
    scale_color_manual(values = my_pallete)+
    scale_fill_manual(values = my_pallete)+
    my_theme+
    theme(legend.position = c(0.2,0.85),
          legend.background = element_blank())+
    labs(x=x_title,color="",fill="")  
  return(g)
}

# This function makes scatter plots of syc metrics vs predictors, color coded by season
# Input include merged syc_df, with predictors merged, in long data form
# the response and predictor variable name
syc_scatter_long <- function(df_long,res_varname,pre_varname,y_title){
  # Calculate R2 for each season
  R2_df <- df_long %>%
    group_by(Season) %>%
    do({
      mod <- lm(.data[[res_varname]] ~ .data[[pre_varname]], data = .)
      s <- summary(mod)
      tibble(
        R2   = s$r.squared,
        pval = coef(s)[2, 4]  # p-value for slope
      )
    }) %>%
    ungroup() %>%
    mutate(
      # format p-values in ggsignif style
      p_label = case_when(
        pval < 0.001 ~ "***",
        pval < 0.01  ~ "**",
        pval < 0.05  ~ "*",
        TRUE         ~ "ns"
      ),
      label = paste0(Season, ": R² = ", round(R2, 2),
                     " (", p_label, ")")
    )
  
  # Get location for labels
  x_max <- max(df_long[[pre_varname]],na.rm=TRUE)
  y_range <- range(df_long[[res_varname]],na.rm=TRUE)
  y_span <- y_range[2] - y_range[1]
  y_top <- y_range[2]
  y_offset <- 0.15 * y_span
  
  # Make scatter plot
  g <- ggplot(data=df_long,aes(x=.data[[pre_varname]],
                               y=.data[[res_varname]],color=Season,fill=Season))+
    geom_point(size=4,shape=21,alpha=0.8,color="black")+
    # Add fitted line
    geom_smooth(method = "lm",se = TRUE,linetype = "solid",alpha = 0.25)+
    scale_color_manual(values = season_color)+
    scale_fill_manual(values = season_color)+
    geom_text(
      data = R2_df,
      aes(x = x_max,
          y = y_top + 0.6 * y_span - (as.numeric(factor(Season))) * y_offset,
          label = label, color = Season),
      hjust = 1, vjust = 1, size = 5, show.legend = FALSE
    ) +
    my_theme+
    labs(x = pre_varname,y = y_title,color="")
  
  return(g)
}

# This function is to bin the proportion of positive values into bins of predictors
# Then make a line plot showing the three seasons
# Input includes the df, with delta_response as in binary (+ or -)
# pre_varname: predictor name
# season_color
prop_sign_plot <- function(df,pre_varname,season_color){
  # Calculate proportion of positive values in each predictor bin
  df <- df %>%
    filter(!is.na(.data[[pre_varname]]),!is.na(delta_response))
  
  # Compute unique quantile breaks
  q_breaks <- unique(quantile(df[[pre_varname]], probs = seq(0, 1, 0.2), na.rm = TRUE))
  
  # If fewer than 2 unique breaks, skip plotting
  if(length(q_breaks) < 2){
    warning(paste("Skipping", pre_varname, "- not enough unique quantile breaks"))
    return(NULL)
  }
  
  df <- df %>%
    mutate(
      # Bin predictors
      pred_bin = cut(.data[[pre_varname]],
                     breaks = q_breaks,
                     include.lowest=TRUE)
    ) %>%
    group_by(pred_bin,Season) %>%
    summarise(Prop_pos = mean(delta_response == "+",na.rm=TRUE),.groups = "drop")
  # Make plot
  g <- ggplot(data=df,aes(x = as.numeric(pred_bin),y=Prop_pos,group=Season,color=Season,fill=Season))+
    geom_line(linewidth = 1)+
    geom_point(shape=21,size=4,color="black")+
    scale_x_continuous(
      breaks = 1:5,
      labels = c("0.2", "0.4", "0.6", "0.8", "1.0"),
      name = paste(pre_varname, "quantiles")
    ) +
    scale_color_manual(values=season_color)+
    scale_fill_manual(values=season_color)+
    labs(x=paste(pre_varname,"quantiles"),
         y="Proportion of + values",
         color="",fill="")+
    my_theme+
    theme(legend.position = c(0.15,0.85),
          legend.background = element_blank())
  return(g)
}

# This function is to make boxplot and violin plots for delta-signs across predictor
prop_sign_box <- function(df,pre_varname,my_color){
  # remove NAs
  df_clean <- df %>% filter(!is.na(.data[[pre_varname]]), !is.na(delta_response))
  
  # run significance test (Wilcoxon) by Season
  sig_df <- df_clean %>%
    group_by(Season) %>%
    summarise(
      p_value = tryCatch(
        wilcox.test(.data[[pre_varname]] ~ delta_response)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(sig_label = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE ~ "ns"
    ))
  
  # find top y position for label placement
  y_pos <- max(df_clean[[pre_varname]], na.rm = TRUE) + 
    (range(df_clean[[pre_varname]])[2] - range(df_clean[[pre_varname]])[1]) * 0.1
  
  g <- ggplot(df_clean,aes(x = delta_response,y=.data[[pre_varname]],
                     fill=delta_response,color=delta_response))+
    geom_half_violin(alpha=0.5,color=NA)+
    geom_boxplot(width=0.1,color="black",outlier.color = NA)+
    geom_jitter(aes(x = as.numeric(as.factor(delta_response))+0.2),
                position = position_jitter(width=0.1),
                alpha=0.7)+
    scale_color_manual(values=my_color)+
    scale_fill_manual(values=my_color)+
    facet_wrap(~Season)+
    geom_text(
      data=sig_df,
      aes(x=1.5,y=y_pos,label=sig_label),
      inherit.aes = FALSE,
      size=5,color="red"
    )+
    my_theme+
    labs(x="",y=pre_varname)
  
  return(g)
}

plot_delta_heatmaps <- function(results_delta_all, fill_name,season,palette_name,direction) {
  g <- ggplot(results_delta_all %>% filter(Season == season),
              aes(x = predictor,y=pair,
                  fill = if_else(pval < 0.05,.data[[fill_name]],NA)))+
    geom_tile(color="white",linewidth = 0.3)+
    scale_fill_distiller(
      palette = palette_name,
      na.value = "grey",
      direction = direction,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black"
      )
    )+
    facet_wrap(~ response,ncol=4)+
    my_theme+
    labs(x = "",y="",fill=fill_name)+
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(season)  
  return(g)
}

# This function is to make chord diagram
# Input include:
# plot_df: df of sychrony metrics between pairs, for both directions. Should have 4 columns
# cols. Colors of the bands
# var_order: order of variables
plot_chord_diagram <- function(plot_df,cols,var_order){
  # Ribbon color = whichever direction is stronger for that pair
  ribbon_cols <- ifelse(
    plot_df$value_direct1 >= plot_df$value_direct2,
    cols[as.character(plot_df$from)],   # color by 'from'
    cols[as.character(plot_df$to)]      # color by 'to'
  )  
  
  # save/restore ONLY what we change (do NOT capture mfrow)
  old_par <- par(c("mar","xaxs","yaxs"))
  on.exit(par(old_par), add = TRUE)
  par(mar = c(2,2,2,2), xaxs = "i", yaxs = "i")
  
  circos.clear()
  circos.par(start.degree = 90, 
             gap.after = rep(4, length(unique(c(plot_df$from,plot_df$to)))),
             track.height = 0.05,
             cell.padding = c(0,0,0,0),
             canvas.xlim = c(-1.2,1.2),
             canvas.ylim = c(-1.2,1.2),
             points.overflow.warning = FALSE)
  
  chordDiagram(
    x = plot_df[,1:4],
    order = var_order,
    grid.col = cols,
    col = ribbon_cols,
    transparency = 0.2,
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
    circos.axis(h = "top", labels = FALSE, 
                major.tick.length = convert_y(1.5, "mm"),
                minor.ticks = 0)
    
    # bigger labels; psi as expression
    lab_map <- list(ET = "ΔET", 
                    VPD = "ΔVPD", 
                    TA = "ΔT", 
                    psi = "Δψ")
    
    circos.text(mean(xlim), ylim[1] + 160, labels = lab_map[[sector]],
                facing = "bending", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 3)
  })
}

# Function to Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# This function is to make correlation matrix plot
# Input is the correlation matrix derived from Hmisc::rcorr
plot_CM <- function(CM,my_label){
  # Get lower r
  CM_r <- get_lower_tri(CM$r)
  # Get lower p
  CM_p <- get_lower_tri(CM$P)
  # Reshape the correlation matrix
  CM_df <- melt(CM_r)
  CM_df$P <- melt(CM_p)$value
  # Remove rows where Var1 = Var2
  CM_df <- CM_df %>%
    filter(Var1 != Var2) %>%
    mutate(P = case_when(
      P < 0.001 ~ "***",
      P < 0.01 ~ "**",
      P < 0.05 ~ "*",
      P >= 0.05 ~ NA
    ),
    r_label = sprintf("%.2f",value),
    label = ifelse(is.na(P),r_label,paste0(r_label,"\n",P))) %>%
    filter(!is.na(value))
  my_color <- rev(RColorBrewer::brewer.pal(11,"RdBu"))
  # Make plot of correlation matrix
  g_CM <- ggplot(CM_df,aes(x=Var1,y=Var2,fill=value))+
    geom_tile(color="white")+
    scale_fill_gradientn(colours = my_color,
                         na.value = "white",
                         limits= c(-1,1))+
    geom_text(aes(label = label))+
    theme(axis.text.x = element_text(angle = 45,hjust=1),
          axis.line=element_line(color="black",size=0.2),
          panel.background = element_blank(),
          text = element_text(size=16),
          axis.text = element_text(size=16),
          panel.border = element_rect(colour="black",fill=NA,size=1),
          legend.frame = element_rect(color="black"),
          legend.ticks = element_line(color="black"),
          legend.title = element_text(margin = margin(b=15)),
          axis.ticks = element_line(size=0.2,color="black"),
          axis.ticks.length = unit(0.05,'in'),
          panel.grid.major = element_line(color="grey"),
          legend.position = c(0.15,0.8),
          legend.background = element_rect(color="black"),
          legend.margin = margin(10,15,15,15))+
    labs(fill = "Spearman ρ",x="",y="")+
    scale_x_discrete(labels=my_label[-1])+
    scale_y_discrete(labels=my_label)
  return(g_CM)
}

# This function is to make box plots
plot_box <- function(df,varname1,varname2,fill_color,x_title,y_title){
  g <-   ggplot(data=df,aes(x=.data[[varname1]],y=.data[[varname2]]))+
    geom_boxplot(fill=fill_color,outlier.color = "grey")+
    my_theme+
    labs(x=x_title,y=y_title)

  return(g)
}

# This function makes boxplot, and show sample size on the side
# df: data frame
# varname1: variable name for the x axis
# varname2: variable name for the y axis
# fill_color: color for the boxes
# show_n: logic. TRUE or FALSE: whether sample number should be shown
# h_just: a negative value for the position of sample size on the right of the panel
# show_letters: logic. TRUE or FALSE: whether letter comparisons should be shown
# letter_offset_frac: letter position
# box_violin decides whether box plot or violin plot should be plotted
plot_box_groups <- function(df,varname1,varname2,fill_color,x_title,y_title,
                            show_n = TRUE,h_just=-4,
                            show_letters = TRUE,letter_offset_frac=0.03,box_violin="Violin"){
  # Keep only complete cases for 
  df <- df[,c(varname1,varname2)]
  df <- df[complete.cases(df),,drop=FALSE]
  # Make sure group is a factor
  #df[[varname2]] <- factor(df[[varname2]])
  # Save factor levels BEFORE subsetting (if varname2 is already a factor)
  lev2 <- if (is.factor(df[[varname2]])) levels(df[[varname2]]) else NULL
  
  # Preserve levels (so blank levels remain)
  if (!is.null(lev2)) {
    df[[varname2]] <- factor(df[[varname2]], levels = lev2)
  } else {
    df[[varname2]] <- factor(df[[varname2]])
  }
  
  
  # Make base plot
  if(box_violin == "Box"){
    g <- ggplot(data=df,aes(x=.data[[varname1]],y=.data[[varname2]]))+
      geom_boxplot(fill=fill_color,outlier.colour = "grey")+
      my_theme+
      labs(x=x_title,y=y_title)+
      # This ensures adding text outside the border
      coord_cartesian(clip = "off")+
      # Make room for texts on the right
      theme(plot.margin = margin(r=80))+
      scale_y_discrete(drop=FALSE)   
  }else{
    g <- ggplot(data=df,aes(x=.data[[varname2]],y=.data[[varname1]]))+
      geom_half_violin(fill=fill_color,alpha=0.5,color=NA,side = "r")+
      geom_boxplot(width=0.1,color="black",outlier.color = NA)+
      geom_jitter(aes(x = as.numeric(.data[[varname2]])-0.2),
                  position = position_jitter(width=0.1),
                  alpha=0.7,
                  color=fill_color)+
      # Flip the coordinate and ensure adding text outside the border
      coord_flip(clip="off")+
      my_theme+
      labs(y=x_title,x=y_title)+
      # Make room for texts on the right
      theme(plot.margin = margin(r=80))+
      scale_x_discrete(drop=FALSE)    
  }

  # Sample size on the right side
  if(show_n){
    n_df <- df %>%
      count(.data[[varname2]],name="n") %>%
      rename(group = !!varname2)
    
    # Add sample size to the plot, outside the border
    if(box_violin == "Box"){
      g <- g +
        geom_text(data=n_df,aes(x=Inf,y=group,label = n),
                  inherit.aes = FALSE,
                  hjust=h_just,
                  size=4)  
    }else{
      g <- g +
        geom_text(data=n_df,aes(y=Inf,x=group,label = n),
                  inherit.aes = FALSE,
                  hjust=h_just,
                  size=4)
    }
  }
  
  # Kruskal wallis test to compare synchrony metrics across groups
  k_result <- kruskal.test(df[[varname1]]~df[[varname2]])
  # Get p-value
  p_txt <- paste0("p = ",formatC(k_result$p.value, format = "e", digits = 2))
  # If not significant, add p-value to the plot
  if(k_result$p.value > 0.05){
    if(box_violin == "Box"){
      g <- g +
        annotate("text",label = p_txt,x=Inf,y=-Inf,hjust=1.1,vjust=-0.8,size=5)      
    }else{
      g <- g +
        annotate("text",label = p_txt,y=Inf,x=-Inf,hjust=1.1,vjust=-0.8,size=5)  
    }
  
  }else {
    # If significant
    # Show significant letters
    if(show_letters){
      # Dunn test
      dunn_res <- FSA::dunnTest(
        x = df[[varname1]],
        g = df[[varname2]],
        method = "bh"
      )$res
      
      # Split comparison column safely
      pairs <- strsplit(dunn_res$Comparison, " - ")
      
      group1 <- sapply(pairs, `[`, 1)
      group2 <- sapply(pairs, `[`, 2)
      
      # Build named p-value vector correctly
      pvec <- dunn_res$P.adj
      names(pvec) <- paste(group1, group2, sep = "-")
      
      # Generate compact letters
      letters <- multcompView::multcompLetters(pvec)$Letters
      
      # check if all letters are the same
      all_same <- length(unique(letters)) == 1
      
      letter_df <- data.frame(
        group = names(letters),
        letters = letters,
        stringsAsFactors = FALSE
      )
      
      # need x_span for letter placement
      x_rng  <- range(df[[varname1]], na.rm = TRUE)
      x_span <- diff(x_rng)
      
      x_max_by_group <- df %>%
        group_by(.data[[varname2]]) %>%
        summarise(x_max = max(.data[[varname1]], na.rm = TRUE), .groups = "drop") %>%
        rename(group = !!varname2)
      
      letter_df <- dplyr::left_join(letter_df, x_max_by_group, by = "group") %>%
        dplyr::mutate(x = x_max + letter_offset_frac * x_span)
      
      # Add letters
      if(all_same){
        # No pairwise differences → show p-value instead
        if(box_violin == "Box"){
          g <- g +
            annotate("text",
                     label = p_txt,
                     x = Inf,
                     y = -Inf,
                     hjust = 1.1,
                     vjust = -0.8,
                     size = 5)
        } else {
          g <- g +
            annotate("text",
                     label = p_txt,
                     y = Inf,
                     x = -Inf,
                     hjust = 1.1,
                     vjust = -0.8,
                     size = 5)
        }
      }else{
        # Add letters (your original code unchanged)
        if(box_violin == "Box"){
          g <- g +
            geom_text(data=letter_df,
                      aes(x=x,y=group,label=letters),
                      inherit.aes = FALSE,
                      hjust=0,
                      size=4)        
        }else{
          g <- g +
            geom_text(data=letter_df,
                      aes(y=x,x=group,label=letters),
                      inherit.aes = FALSE,
                      vjust=0,
                      size=4)
        }
      }
    }  
  }
  
  return(g)
}

# This function is to make Peak TE vs lag plots
plot_TE_vs_lag <- function(shape_var){
  g <- ggplot(data=syc_df,aes(x = GS_best_lag,y=forcats::fct_rev(source_sink)))+
    geom_point(aes(size=GS_daily_p_TE^4,color=source_sink,shape=.data[[shape_var]]))+
    scale_color_manual(values = pair_color)+
    my_theme+
    labs(x = "Lag (h)",y="")+
    scale_y_discrete(labels = c(
      psi_ET = expression(psi %->% ET),
      ET_psi = expression(ET %->% psi),
      VPD_ET = expression(VPD %->% ET),
      ET_VPD = expression(ET %->% VPD),
      TA_ET  = expression(T[air] %->% ET),
      ET_TA  = expression(ET %->% T[air])
    ))+
    scale_size_continuous(
      name = "Peak TE (%)",
      breaks = c(2, 4, 8, 16)^4,   # breaks in plotted (transformed) space
      labels = c("2", "4", "8", "16"),
      range = c(2, 10)
    )+
    guides(color = "none")+
    theme(legend.position = "right")
  return(g)
}

# Vertical version of box plot distributions across groups, specifically designed for lags
plot_lag_violin_box <- function(df, group_var, fill_color = "grey80",
                                show_n = TRUE,
                                show_letters = TRUE,
                                letter_offset_frac = 0.03) {
  
  # keep needed columns
  df <- df[, c("GS_best_lag", "source_sink", group_var)]
  df <- df[complete.cases(df), , drop = FALSE]
  
  # preserve group levels if already factor
  lev_group <- if (is.factor(df[[group_var]])) levels(df[[group_var]]) else NULL
  if (!is.null(lev_group)) {
    df[[group_var]] <- factor(df[[group_var]], levels = lev_group)
  } else {
    df[[group_var]] <- factor(df[[group_var]])
  }
  
  # preserve source_sink levels if already factor
  lev_ss <- if (is.factor(df$source_sink)) levels(df$source_sink) else NULL
  if (!is.null(lev_ss)) {
    df$source_sink <- factor(df$source_sink, levels = lev_ss)
  } else {
    df$source_sink <- factor(df$source_sink)
  }
  
  # facet labels as expressions
  facet_labs <- c(
    psi_ET = "psi %->% ET",
    ET_psi = "ET %->% psi",
    VPD_ET = "VPD %->% ET",
    ET_VPD = "ET %->% VPD",
    TA_ET  = "T[air] %->% ET",
    ET_TA  = "ET %->% T[air]"
  )
  
  # base plot
  g <- ggplot(df, aes(x = .data[[group_var]], y = GS_best_lag)) +
    geom_violin(fill = fill_color, alpha = 0.5, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.color = NA) +
    facet_wrap(
      ~source_sink,
      scales = "free_y",
      labeller = as_labeller(facet_labs, label_parsed)
    ) +
    labs(x = NULL, y = "Lag (h)") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    coord_cartesian(clip = "off") +
    my_theme
  
  # sample size
  if (show_n) {
    n_df <- df %>%
      dplyr::group_by(source_sink, .data[[group_var]]) %>%
      dplyr::summarise(
        n = dplyr::n(),
        y_max = max(GS_best_lag, na.rm = TRUE),
        y_min = min(GS_best_lag, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        y = y_max + 0.08 * ifelse(y_max - y_min == 0, 1, y_max - y_min)
      )
    
    g <- g +
      geom_text(
        data = n_df,
        aes(x = .data[[group_var]], y = y, label = n),
        inherit.aes = FALSE,
        size = 4
      )
  }
  
  # store annotations separately
  p_list <- list()
  letter_list <- list()
  
  for (ss in levels(df$source_sink)) {
    sub_df <- df[df$source_sink == ss, , drop = FALSE]
    
    # remove unused factor levels within facet
    sub_df[[group_var]] <- droplevels(sub_df[[group_var]])
    
    # skip if too few groups
    if (nrow(sub_df) == 0 || length(unique(na.omit(sub_df[[group_var]]))) < 2) next
    
    # Kruskal-Wallis
    k_result <- kruskal.test(sub_df$GS_best_lag ~ sub_df[[group_var]])
    p_txt <- paste0("p = ", formatC(k_result$p.value, format = "e", digits = 2))
    
    y_rng  <- range(sub_df$GS_best_lag, na.rm = TRUE)
    y_span <- diff(y_rng)
    if (y_span == 0) y_span <- 1
    
    first_group <- levels(sub_df[[group_var]])[1]
    
    if (k_result$p.value > 0.05) {
      p_list[[length(p_list) + 1]] <- data.frame(
        source_sink = ss,
        x = first_group,
        y = y_rng[2] + 0.15 * y_span,
        label = p_txt,
        stringsAsFactors = FALSE
      )
      
    } else if (show_letters) {
      dunn_res <- FSA::dunnTest(
        x = sub_df$GS_best_lag,
        g = sub_df[[group_var]],
        method = "bh"
      )$res
      
      pairs <- strsplit(dunn_res$Comparison, " - ")
      group1 <- sapply(pairs, `[`, 1)
      group2 <- sapply(pairs, `[`, 2)
      
      pvec <- dunn_res$P.adj
      names(pvec) <- paste(group1, group2, sep = "-")
      
      letters <- multcompView::multcompLetters(pvec)$Letters
      all_same <- length(unique(letters)) == 1
      
      if (all_same) {
        p_list[[length(p_list) + 1]] <- data.frame(
          source_sink = ss,
          x = first_group,
          y = y_rng[2] + 0.15 * y_span,
          label = p_txt,
          stringsAsFactors = FALSE
        )
      } else {
        letter_df <- data.frame(
          group = names(letters),
          label = letters,
          stringsAsFactors = FALSE
        )
        
        y_max_df <- sub_df %>%
          dplyr::group_by(.data[[group_var]]) %>%
          dplyr::summarise(
            y_max = max(GS_best_lag, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          dplyr::rename(group = !!group_var)
        
        letter_df <- dplyr::left_join(letter_df, y_max_df, by = "group") %>%
          dplyr::mutate(
            source_sink = ss,
            x = group,
            y = y_max + letter_offset_frac * y_span
          )
        
        letter_list[[length(letter_list) + 1]] <- letter_df
      }
    }
  }
  
  if (length(p_list) > 0) {
    p_df <- dplyr::bind_rows(p_list)
    g <- g +
      geom_text(
        data = p_df,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        size = 4
      )
  }
  
  if (length(letter_list) > 0) {
    letter_df <- dplyr::bind_rows(letter_list)
    g <- g +
      geom_text(
        data = letter_df,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        vjust = 0,
        size = 4
      )
  }
  
  return(g)
}