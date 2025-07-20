# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# This code is to make plots

library(ggplot2)
library(cowplot)
library(RColorBrewer)

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
# varcolor: color for this line
TS_all <- function(varname,df,y_title,varcolor){
  g <- ggplot(data=df,aes(x = Time,y=.data[[varname]]))+
    geom_line(color=varcolor)+
    my_theme+
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
        mean = mean(.data[[varname]],na.rm=TRUE))
  }else if(cycle == "Diurnal"){
    # For diurnal cycle
    df_tmp <- df %>%
      mutate(Hour = hour(Time)) %>%
      # Calculate hourly mean across the day
      group_by(Hour) %>%
      summarise(
        Time = first(Hour),
        mean = mean(.data[[varname]],na.rm=TRUE))
  }
  return(df_tmp)
}

# This function is to make annual or diurnal time series TS plots
# Input includes:
# varname: the variable name in the df
# df_cycle: summarized df of annual or diurnal cycle
# cycle: "Annual" or "Diurnal"
TS_cycle <- function(df_cycle,cycle){
  df_cycle$Type <- factor(df_cycle$Type,levels=c("Original","Diurnal mean","Diurnal anomaly"))
  g <- ggplot(df_cycle,aes(x=Time,y=mean,color=Type))+
    geom_line(aes(y=mean),linewidth=1,alpha=0.7)+
    my_theme+
    scale_color_manual(values = my_color)+
    labs(x = "",y=y_title,color="")
  if(cycle == "Annual"){
    g <- g+
      scale_x_date(date_breaks = "2 month",date_labels = "%b")
  }
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
  g_annual <- TS_cycle(df_annual,"Annual")
  
  # Diurnal cycle
  # Make a df for diurnal cycle data
  df_diurnal <- rbind(var_cycle(varname,df,"Diurnal"),
                      var_cycle(varname_mean,df,"Diurnal"),
                      var_cycle(varname_anomaly,df,"Diurnal"))
  df_diurnal$Type <- rep(c("Original",
                           "Diurnal mean",
                           "Diurnal anomaly"),
                         each = nrow(df_diurnal)/3)
  g_diurnal <- TS_cycle(df_diurnal,"Diurnal")
  
  # Put all figures together
  g1 <- plot_grid(plotlist = g_hist_original_ls,nrow=1)
  g2 <- plot_grid(plotlist = g_hist_mean_ls,nrow=1)
  g3 <- plot_grid(plotlist = g_hist_anomaly_ls,nrow=1)
  g4 <- plot_grid(g_annual,g_diurnal,nrow=1)
  g_all <- plot_grid(g_original,g1,
                  g_mean,g2,
                  g_anomaly,g3,
                  ncol=1,align="v",axis = "lr")
  g_all <- plot_grid(g_all,g4,nrow=2,rel_heights = c(6,1),align="hv")
  return(g_all)
}

# This function is to plot TE vs lag from var1 to var2
# Input is the TE_df
# And the name of the two variables varname1 and varname2
# And title of the plot
# norm: "None","Shannon" or "Theory"
TE_lag_plot <- function(TE_df,title,norm){
  if(norm == "None"){
    y_title <- "TE (bits)"
  }else{
    y_title <- "TE (%)"
  }
  # Plot of TE
  g_TE <- ggplot(data=TE_df,aes(x=Lag,y=TE))+
    geom_line(linewidth = 0.4)+
    #geom_ribbon(aes(ymin = TE-TE_se,ymax=TE+TE_se),alpha=0.3,color=NA)+
    geom_line(aes(y=TEcrit),linewidth=0.4,linetype = "dashed",color="royalblue2")+
    my_theme+
    labs(x = "Lag (days)",y=y_title,color="",fill="")+
    ggtitle(title)
  return(g_TE)
}



