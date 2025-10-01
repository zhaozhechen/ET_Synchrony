# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# This code is to make plots

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(gghalves)
library(stringr)
library(sf)

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
  legend.position = "none"
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
syc_map <- function(df,varname,palette_name,legend_title,g_title,color_limits){
  g <- ggplot()+
    geom_sf(data=CONUS,fill="grey",color="black",alpha=0.3)+
    geom_point(data=df,aes(x=longitude,y=latitude,
                           fill=.data[[varname]]),
               size=5,alpha=0.8,shape=21,color="black")+
    scale_fill_distiller(
      palette = palette_name, direction = -1,
      limits = color_limits,
      oob = scales::squish,  # keeps values inside limits
      breaks = seq(color_limits[1], color_limits[2], length.out = 5),
      labels = function(x){
        labs <- as.character(x)
        labs[1] <- paste0("≤ ", labs[1])   # add ≤ to min
        labs[length(labs)] <- paste0("≥ ", labs[length(labs)]) # add ≥ to max
        labs
      },
      guide = guide_colorbar(
        ticks = TRUE,
        frame.colour = "black",
        ticks.colour = "black"
      )
    )+
    labs(fill = legend_title)+
    ggtitle(g_title)+
    map_theme
  return(g)
}

# This is the same map function but for discrete values
syc_map_disc <- function(df,varname,legend_title,g_title){
  g <- ggplot()+
    geom_sf(data=CONUS,fill="grey",color="black",alpha=0.3)+
    geom_point(data=df,aes(x=longitude,y=latitude,
                           fill=.data[[varname]]),
               size=5,alpha=0.8,shape=21,color="black")+
    scale_fill_manual(
      values = my_color,
      guide = guide_legend(override.aes = list(size = 4))
    )+
    labs(fill = legend_title)+
    ggtitle(g_title)+
    map_theme+
    theme(legend.position = "bottom")
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
  p <- signif(summary(lm)$coefficients[2,4],2)
  R2 <- round(summary(lm)$r.squared,3)
  
  g <- ggplot(data=df,aes(x=.data[[varname1]],y=.data[[varname2]],color=.data[[group_name]]))+
    geom_point(size=2,alpha=0.7)+
    my_theme+
    labs(x=xtitle,y=ytitle,color="")+
    geom_abline(slope = k,intercept = b)+
    annotate(
      geom = "text",
      x=Inf,y=-Inf,
      label = bquote(italic(p) == .(p)),
      hjust=1.2,vjust=-2,
      size=5)+
    annotate(
      geom = "text",
      x=Inf,y=-Inf,
      label = bquote(R^2 == .(R2)),
      hjust=1.2,vjust=-4,
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
