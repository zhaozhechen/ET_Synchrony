# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# These codes are to analyze synchrony metrics based on TE results

# This function calculates synchrony metrics
# Input include:
# TE_df: full TE_df from the source to the sink
# Output Synchrony metrics including
# p_TE: peak daily TE
# p_lag: lag corresponding to peak TE
# agg_TE: daily aggregated TE in the first 24 hours
# memory: four options (mem1,mem3,mem4,and,mem5). mem2 is not calculated
cal_syc_metrics <- function(TE_df){
  # Normalize TE and TEcrit
  TE_df_tmp <- TE_df %>%
    mutate(TE_norm = TE/Hy * 100,
           TEcrit_norm = TEcrit/Hy * 100,
           # Add significance flag
           sig = if_else(!is.na(TE_norm) & !is.na(TEcrit_norm) & (TE_norm > TEcrit_norm),
                         TRUE,FALSE)) %>%
    select(Lag,TE_norm,TEcrit_norm,sig)
  
  # Extract synchrony metrics --------
  # Only consider the first 24 hours
  TE_df_24h <- TE_df_tmp[1:24,]
  # if all TE values are insignificant
  if(sum(TE_df_24h$sig,na.rm=TRUE)==0){
    # All synchrony metrics are NA
    syc_metrics <- rep(NA,8)
  }else{
    # Daily peak TE
    p_TE <- max(TE_df_24h$TE_norm[TE_df_24h$sig],na.rm=TRUE)
    # Get corresponding lag
    p_lag <- TE_df_24h$Lag[which(TE_df_24h$TE_norm == p_TE)[1]]
    # Calculate aggregate TE within the first 24 hours
    agg_TE <- sum(TE_df_24h$TE_norm[TE_df_24h$sig],na.rm=TRUE)
    
    # Calculate memory ---------
    # Option 1: From lag 0 to first non-significant lag
    nonsig_lag1 <- TE_df_tmp %>% filter(sig == FALSE)
    mem1 <- if(nrow(nonsig_lag1) > 0) nonsig_lag1$Lag[1] else max(TE_df_tmp$Lag)
    
    # Option 2: Width of significant TE that includes peak TE
    peak_idx <- which(TE_df_tmp$Lag == p_lag)
    # Walk left from peak
    left_idx <- peak_idx
    while(left_idx > 1 && TE_df_tmp$sig[left_idx - 1]) {
      left_idx <- left_idx - 1
    }
    # Walk right from peak
    right_idx <- peak_idx
    while(right_idx < nrow(TE_df_tmp) && TE_df_tmp$sig[right_idx + 1]) {
      right_idx <- right_idx + 1
    }
    # Width
    mem2 <- TE_df_tmp$Lag[right_idx] - TE_df_tmp$Lag[left_idx] + 1
    
    # Option 3: From peak to first non-significant lag after the peak
    after_peak <- TE_df_tmp %>% filter(Lag > p_lag)
    nonsig_after_peak <- after_peak$Lag[which(after_peak$sig == FALSE)]
    mem3 <- if(length(nonsig_after_peak) > 0) nonsig_after_peak[1] - p_lag else max(TE_df_tmp$Lag)
    
    # Option 4: Total significant duration (cumulative hours)
    mem4 <- sum(TE_df_tmp$sig,na.rm=TRUE)-1
    
    # Option 5: From lag 0 to first lag after the peak
    mem5 <- if(length(nonsig_after_peak) > 0) nonsig_after_peak[1] else max(TE_df_tmp$Lag)

    # Combine all synchrony metrics -----------
    syc_metrics <- c(p_TE,p_lag,agg_TE,mem1,mem2,mem3,mem4,mem5)
  }
  names(syc_metrics) <- c("pTE","plag","aggTE","mem1","mem2","mem3","mem4","mem5")
  return(syc_metrics)
}

# This function calculates syc metrics for all pairs of variables
# For the selected time period at this site
# Input include:
# Site_ID: Site_ID for the site
# type: should be a character: "full_TS","GS", or "NGS"
# Output: a vector of length 84, for syc metrics for all pairs of variables
cal_syc_metrics_all_pairs <- function(Site_ID,type,var_comb){
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
      syc_metrics <- rep(NA,8)
      names(syc_metrics) <- c("pTE","plag","aggTE","mem1","mem2","mem3","mem4","mem5")
    }else{
      # Calculate Synchrony metrics if valid
      syc_metrics <- cal_syc_metrics(TE_df)
    }
    # Name these metrics
    names(syc_metrics) <- paste0(names(syc_metrics),"_",TE_df_name)
    # Combine all syc metrics
    syc_metrics_all <- c(syc_metrics_all,syc_metrics)
  }
  return(syc_metrics_all)
}

# This function extracts synchrony metrics for each site, 
# for full_TS, GS, and NGS
# Output a df of three rows
syc_site <- function(Site_ID,var_comb){
  # For full year
  syc_metrics_full_TS <- cal_syc_metrics_all_pairs(Site_ID,"full_TS",var_comb)
  # For GS
  syc_metrics_GS <- cal_syc_metrics_all_pairs(Site_ID,"GS",var_comb)
  # For NGS
  syc_metrics_NGS <- cal_syc_metrics_all_pairs(Site_ID,"NGS",var_comb)
  
  # Add Site_ID and GS
  syc_metrics_site <- data.frame(Site_ID = rep(Site_ID,3),
                                 GS = c("FT","GS","NGS"),
                                 rbind(syc_metrics_full_TS,
                                       syc_metrics_GS,
                                       syc_metrics_NGS))
  
  rownames(syc_metrics_site) <- NULL
  return(syc_metrics_site)
}





