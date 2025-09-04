# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# These codes are to analyze synchrony metrics based on TE results

library(lomb) # This is for Lomb-Scargle Periodogram


# This function detects the frequency of TE vs lag plot
# Input: 
# TE_df_tmp: include two columns, Time and normalized TE
# This function outputs the significant period, and its LS plot
# If no significant period, returns NA
get_period <- function(TE_df_tmp){
  TE_df_tmp <- TE_df_tmp %>%
    select(Lag,TE_norm)
  # Use Lomb-Scargle Periodogram to detect
  results <- lsp(TE_df_tmp,from =1/48,to=1/4,ofac =3,type="frequency",alpha=0.01)

  # Make the Lomb-Scargle Periodogram
  g_LSP <- plot_LSP(results)
  
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
# Output include:
# p_TE: peak TE
# p_lag: lag corresponding to peak TE
# memory: time required for TE to drop below critical TE
cal_syc_metrics <- function(TE_df){
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
  LS_results <- get_period(TE_df_tmp)

  # Extract synchrony metrics --------
  # Only consider the first 24 hours
  TE_df_24h <- TE_df_tmp[1:24,]
  # if all TE values are insignificant
  if(sum(TE_df_24h$sig,na.rm=TRUE)==0){
    p_TE <- NA
    p_lag <- NA
    mem <- NA
    agg_TE <- NA
  }else{
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
  }
  
  # Output synchrony metrics
  syc_metrics <- c(p_TE,p_lag,mem,agg_TE,period)
  
  
  
  
  
  # Should directly output LS plot, together with TE plot
  # File name should be passed
  
  
  out <- list(LS_g = LS_g,
              syc_metrics = c(p_TE,p_lag,mem,agg_TE,period))
  return(out)
}

# This function calculates syc metrics for all pairs of variables
# For the selected time period at this site
# Input include:
# Site_ID: Site_ID for the site
# type: should be a character: "full_TS","GS", or "NGS"
# Output: a vector of length 36, for syc metrics for all pairs of variables
Site_ID <- Site_ID
type <- "full_TS"
cal_syc_metrics_all_pairs <- function(Site_ID,type){
  # Get file name for the TE_df_ls
  file_name <- paste0("TE_df_ls_",type,"_",Site_ID,".rds")
  # Read in TE_df_ls
  TE_df_ls <- readRDS(paste0(TE_df_path,file_name))
  # Loop over each variable pair
  # Initialize a vector to store all output syc metrics for all pairs
  syc_metrics_all <- c()
  for(i in 1:nrow(var_comb)){
    # Source and sink name
    source_name <- as.character(var_comb$from[i])
    sink_name <- as.character(var_comb$to[i])
    # Get TE_df for this pair
    TE_df_name <- paste0(source_name,"_to_",sink_name)
    TE_df <- TE_df_ls[[TE_df_name]]
    # Check if this TE_df is valid
    if(nrow(TE_df) < max_lag){
      # All NA
      syc_metrics <- c(NA,NA,NA,NA,NA)
    }else{
      # Calculate Synchrony metrics if valid
      syc_results <- cal_syc_metrics(TE_df,TE_df_name,type)
      
      
      
      
      
      
      
      syc_metrics <- cal_syc_metrics(TE_df)[2]
    }
    # Name these metrics
    names(syc_metrics) <- paste0(c("p_TE_","p_lag_","mem_","agg_TE_"),TE_df_name)
    syc_metrics_all <- c(syc_metrics_all,syc_metrics)
  }
  return(syc_metrics_all)
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
