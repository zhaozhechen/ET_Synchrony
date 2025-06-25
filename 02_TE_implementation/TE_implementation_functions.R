# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# This code includes functions for TE implementation, and data processing before TE implementation

# This function calculates the difference between each two time steps in a time series
# Input df: the data frame of time series variables
# Input varname: variable name to process
delta_TS <- function(df,varname){
  # TS at t0
  TS0 <- head(df[[varname]],-1)
  # TS with 1 time step shifted
  TS1 <- tail(df[[varname]],-1)
  # Get the difference
  delta_var <- TS1 - TS0
  return(delta_var)
}

# This function calculates TE from var1 to var2, for lags from 0 to max_lag
# There are two ways of normalization: Shannon (default) or Theory (which uses logn_bin)
# Input includes:
# var1: source variable
# var2: sink variable
# max_lag: The maximum lag for running TE
# norm: "Shannon" or "Theory"
# alpha: the alpha level used, default 0.05
# Number of bootstrapping, default 300
Cal_TE <- function(var1,var2,max_lag,norm = "Shannon",alpha=0.05,nboot = 300){
  # Initialize a progress bar
  pb <- progress_bar$new(
    format = "Calculating TE[:bar]:percent | Lag :current/:total",
    total = max_lag + 1,clear = FALSE,width=60
  )
  
  # Initialize a list to store all TE
  TE_list <- lapply(0:max_lag, function(lag){
    # Align TS
    var2_lagged <- var2[(lag+1):length(var2)]
    var1_aligned <- var1[1:(length(var1)-lag)]
    
    # Calculate TE(var1,var2),unit: bits
    TE_results <- transfer_entropy(x = var1_aligned,
                                   y = var2_lagged,
                                   lx=1,ly=1, # Markov order (memory of the two variables)
                                   shuffles = 100, # Shuffle to calculate effective TE (for bias correction)
                                   type = "bins", # Method used for discretization of continuous data
                                   bins = n_bin, # Number of bins for discretization
                                   nboot = nboot, # Bootstrapping 300 times, for statistical inference
                                   seed = my_seed, # Set a seed for random number generator
                                   quiet = TRUE)
    # Extract required values
    # TE
    TE <- TE_results$coef[1,1]
    # Effective TE
    ETE <- TE_results$coef[1,2]
    # Standard error for TE (from bootstrapping)
    TE_se <- TE_results$coef[1,3]
    # Quantile of bootstrapping results
    # Using 95th quantile as the critical TE values, from X -> Y
    cr_TE <- quantile(TE_results$boot[1,],1-alpha)
    
    # Calculate values for normalization
    if(norm == "Shannon"){
      var2_disc <- discretize(na.omit(var2_lagged),numBins = n_bin)
      H_value <- entropy(var2_disc,unit="log2")
    }else{
      H_value <- log(n_bin)
    }
    
    # Normlize TE
    TE_df <- data.frame(Lag = lag,
                        TE = TE/H_value*100,
                        ETE = ETE/H_value*100,
                        TE_se = TE_se/H_value*100,
                        cr_TE = cr_TE/H_value*100)
    #message(paste("Complete Lag",lag))
    # Update progress bar
    pb$tick()
    return(TE_df)
  })
  # Combine df for each lag
  TE_df <- bind_rows(TE_list)
  return(TE_df)
}

