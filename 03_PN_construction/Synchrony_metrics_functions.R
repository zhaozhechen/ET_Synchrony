# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# These codes are to analyze synchrony metrics based on TE results

# This function calculates synchrony metrics
# Input include the TE_df
# Output include:
# p_TE: peak TE
# p_lag: lag corresponding to peak TE
# memory: time required for TE to drop below critical TE
cal_syc_metrics <- function(TE_df){
  # Normalize TE and TEcrit
  TE_df_tmp <- TE_df %>%
    mutate(TE_norm = TE/Hy * 100,
           TEcrit_norm = TEcrit/Hy * 100) %>%
    select(Lag,TE_norm,TEcrit_norm)
  # Only consider the first 24 hours
  TE_df_1day <- TE_df_tmp[1:24,]
  # Only consider significant values
  TE_df_1day <- TE_df_1day %>%
    mutate(TE_norm = if_else(TE_norm < TEcrit_norm,
                             NA,TE_norm))
  # Get peak TE
  # if all TE values are insignificant
  if(sum(TE_df_1day$TE_norm,na.rm=TRUE)==0){
    p_TE <- NA
    p_lag <- NA
    memory <- NA
  }else{
    p_TE <- max(TE_df_1day$TE_norm,na.rm=TRUE)
    # Get corresponding lag
    p_lag <- TE_df_1day$Lag[which(TE_df_1day$TE_norm == p_TE)]
    # Get memory
    memory <- TE_df_tmp$Lag[which(TE_df_tmp$TE_norm < TE_df_tmp$TEcrit_norm)]
    # memory has to be after peak TE
    memory <- memory[memory > p_lag][1]
    
    # If TE never drops below critical TE, assign max_lag to it
    if(is.na(memory)){
      memory <- max_lag
    }    
  }
  
  return(c(p_TE,p_lag,memory))
}

# This function calculates syc metrics for all pairs of variables
# For the selected time period at this site
# Input include:
# file_name: This should be the name for the TE_df_ls
# Output: a vector of length 36, for syc metrics for all pairs of variables
cal_syc_metrics_all_pairs <- function(file_name){
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
      syc_metrics <- c(NA,NA,NA)
    }else{
      # Calculate Synchrony metrics if valid
      syc_metrics <- cal_syc_metrics(TE_df)
    }
    # Name these metrics
    names(syc_metrics) <- paste0(c("p_TE_","p_lag_","mem_"),TE_df_name)
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
