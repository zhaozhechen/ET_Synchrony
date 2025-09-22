# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# These codes are to analyze synchrony metrics based on TE results



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
    syc_metrics <- rep(NA,7)
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
    
    # Option 2: Width of TE peak (only significant). Not calculated
    
    # Option 3: From peak to first non-significant lag after the peak
    after_peak <- TE_df_tmp %>% filter(Lag > p_lag)
    nonsig_after_peak <- after_peak$Lag[which(after_peak$sig == FALSE)]
    mem3 <- if(length(nonsig_after_peak) > 0) nonsig_after_peak[1] - p_lag else max(TE_df_tmp$Lag)
    
    # Option 4: Total significant duration (cumulative hours)
    mem4 <- sum(TE_df_tmp$sig,na.rm=TRUE)
    
    # Option 5: From lag 0 to first lag after the peak
    mem5 <- if(length(nonsig_after_peak) > 0) nonsig_after_peak[1] else max(TE_df_tmp$Lag)

    # Combine all synchrony metrics -----------
    syc_metrics <- c(p_TE,p_lag,agg_TE,mem1,mem3,mem4,mem5)
  }
  names(syc_metrics) <- c("p_TE","p_lag","agg_TE","mem1","mem3","mem4","mem5")
  return(syc_metrics)
}


