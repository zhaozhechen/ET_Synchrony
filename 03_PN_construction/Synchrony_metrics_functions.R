# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# These codes are to analyze synchrony metrics based on TE results

library(lomb) # This is for Lomb-Scargle Periodogram


# This function detects the frequency of TE vs lag plot
# Input: 
# TE_df_tmp: include two columns, Time and normalized TE
# This function outputs the significant period, and its LS plot
# If no significant period, returns NA
get_period <- function(TE_df_tmp,p_color){
  TE_df_tmp <- TE_df_tmp %>%
    select(Lag,TE_norm)
  # Use Lomb-Scargle Periodogram to detect
  results <- lsp(TE_df_tmp,from =1/48,to=1/4,ofac =3,type="frequency",alpha=0.05)

  # Make the Lomb-Scargle Periodogram
  g_LSP <- plot_LSP(results,p_color)
  
  # Check if there is a statistically significant peak
  if(!is.null(results$peak.at) && (results$peak > results$sig.level)){
    period <- 1/results$peak.at[1]
  }else{
    period <- NA
  }
  return(list(g = g_LSP,period=period))
}

# This function calculates synchrony metrics
# Input include:
# TE_df: full TE_df from the source to the sink
# m_color: a color vector of five for the plots
# TE_g_title: title for the TE vs lag plot
# Output include a list of two elements:
# g: the combined plots of TE vs lag and its corresponding Lomb Scargle Periodogram
# Synchrony metrics including
# p_TE: peak daily TE
# p_lag: lag corresponding to peak TE
# memory: time required for TE to drop from peak to below critical TE
# agg_TE: daily aggregated TE in the first 24 hours
cal_syc_metrics <- function(TE_df,m_color,TE_g_title){
  # Normalize TE and TEcrit
  TE_df_tmp <- TE_df %>%
    mutate(TE_norm = TE/Hy * 100,
           TEcrit_norm = TEcrit/Hy * 100,
           # Add significance flag
           sig = if_else(!is.na(TE_norm) & !is.na(TEcrit_norm) & (TE_norm > TEcrit_norm),
                         TRUE,FALSE)) %>%
    select(Lag,TE_norm,TEcrit_norm,sig)
  
  # Detect periodicity, no need to consider significance of TE values
  # Use Lomb-Scargle Periodogram to detect
  LS_results <- get_period(TE_df_tmp,m_color[4])

  # Extract synchrony metrics --------
  # Only consider the first 24 hours
  TE_df_24h <- TE_df_tmp[1:24,]
  # if all TE values are insignificant
  if(sum(TE_df_24h$sig,na.rm=TRUE)==0){
    p_TE <- NA
    p_lag <- NA
    mem <- NA
    agg_TE <- NA
    # Also make TE vs lag plot, but no metrics annotated
    g_TE <- ggplot(TE_df_tmp,aes(x=Lag,y=TE_norm))+
      geom_line(linewidth=0.8)+
      labs(x="Lag (h)",y="Uncertainty reduction (%)")+
      geom_line(aes(y=TEcrit_norm),
                linewidth = 0.8,
                linetype = "dashed",
                color=m_color[2])+
      my_theme
  }else{
    # Daily maximum TE value
    p_TE <- max(TE_df_24h$TE_norm[TE_df_24h$sig],na.rm=TRUE)
    # Get corresponding lag
    p_lag <- TE_df_24h$Lag[which(TE_df_24h$TE_norm == p_TE)[1]]
    
    # Find the first non-significant lag after the peak
    TE_after_peak <- TE_df_tmp %>%
      filter(Lag > p_lag)
    nonsig_lags <- TE_after_peak$Lag[which(TE_after_peak$sig == FALSE)]
    
    # Get memory
    if(length(nonsig_lags) == 0){
      mem <- max_lag
    }else{
      mem <- nonsig_lags[1] - p_lag
    }
    
    # Calculate aggregate TE within the first 24 hours
    agg_TE <- sum(TE_df_24h$TE_norm[TE_df_24h$sig],na.rm=TRUE)
    
    # Make the TE_lag plot, with synchorny metrics labeled
    g_TE <- TE_norm_lag_plot(TE_df_tmp,p_lag,p_TE,mem,m_color)
  }
  
  # synchrony metrics for output
  syc_metrics <- c(p_TE,p_lag,mem,agg_TE,LS_results$period)
  # Combine TE plot and LS plot
  g_TE <- g_TE +
    ggtitle(TE_g_title)
  g <- plot_grid(g_TE,LS_results$g,
                 nrow=1,align = "hv",
                 rel_widths = c(1.5,1))

  out <- list(g,syc_metrics)
  return(out)
}

# This function calculates syc metrics for all pairs of variables
# For the selected time period at this site
# Input include:
# Site_ID: Site_ID for the site
# type: should be a character: "full_TS","GS", or "NGS"
# m_color: a color of five for plotting
# Output: a vector of length 36, for syc metrics for all pairs of variables
cal_syc_metrics_all_pairs <- function(Site_ID,type,m_color){
  # Get file name for the TE_df_ls
  file_name <- paste0("TE_df_ls_",type,"_",Site_ID,".rds")
  # Read in TE_df_ls
  TE_df_ls <- readRDS(paste0(TE_df_path,file_name))
  # Loop over each variable pair
  # Initialize a vector to store all output syc metrics for all pairs
  syc_metrics_all <- c()
  # Initialize a list to store plots for all pairs
  g_all <- list()
  for(i in 1:nrow(var_comb)){
    # Source and sink name
    source_name <- as.character(var_comb$from[i])
    sink_name <- as.character(var_comb$to[i])
    # Get TE_df for this pair
    TE_df_name <- paste0(source_name,"_to_",sink_name)
    TE_df <- TE_df_ls[[TE_df_name]]
    # Get title for the plot
    TE_g_title <- bquote(Delta~.(as.name(source_name))~"\u2192"~Delta~.(as.name(sink_name)))
    # Check if this TE_df is valid
    if(nrow(TE_df) < max_lag){
      # All NA
      syc_metrics <- c(NA,NA,NA,NA,NA)
      # Make a blank plot as place holder
      g <- ggplot()+
        xlim(0,5)+
        ylim(0,5)+
        ggtitle(TE_g_title)+
        my_theme
    }else{
      # Calculate Synchrony metrics if valid
      syc_results <- cal_syc_metrics(TE_df,m_color,TE_g_title)
      syc_metrics <- syc_results[[2]]
      g <- syc_results[[1]]
    }
    # Name these metrics
    names(syc_metrics) <- paste0(c("p_TE_","p_lag_","mem_","agg_TE_","period_"),TE_df_name)
    syc_metrics_all <- c(syc_metrics_all,syc_metrics)
    # Add the figure to the list
    g_all[[i]] <- g
  }
  # Combine all plots and output
  g_out <- plot_grid(plotlist = g_all,
                 ncol=3,
                 align="hv")
  # Output this figure
  print_g(g_out,paste0("TE_LS_plots_",type,"_",Site_ID),
          16,12)
  
  return(syc_metrics_all)
}

# This function get synchrony metrics for each site, 
# and output TE vs lag + LSP plots for each site
syc_site <- function(Site_ID){
  # For full year
  syc_metrics_full_TS <- cal_syc_metrics_all_pairs(Site_ID,"full_TS",my_color)
  # For GS
  syc_metrics_GS <- cal_syc_metrics_all_pairs(Site_ID,"GS",my_color)
  # For NGS
  syc_metrics_NGS <- cal_syc_metrics_all_pairs(Site_ID,"NGS",my_color)
  
  # Combine these metrics together
  syc_metrics_site <- as.data.frame(rbind(syc_metrics_full_TS,
                                          syc_metrics_GS,
                                          syc_metrics_NGS))
  # Add Site ID
  syc_metrics_site$Site_ID <- Site_ID
  # Move Site_ID to the first column
  syc_metrics_site <- syc_metrics_site[, c("Site_ID", setdiff(names(syc_metrics_site), "Site_ID"))]
  # Add time period: full-TS, GS, or NGS
  syc_metrics_site$GS <- c("FT","GS","NGS")
  rownames(syc_metrics_site) <- NULL
  return(syc_metrics_site)
}

# This function is to conduct Kruskal-Wallis test to compare syc metrics across group
# df
# y_varname: variable name to test
# group_name: variable name that the variables need to be grouped by
syc_compare <- function(df,y_varname,group_name){
  f <- as.formula(paste0(y_varname,"~",group_name))
  test <- kruskal.test(f,data=df)
  p_value <- signif(test$p.value,2)
  return(p_value)
}
