# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# This code includes core functions for TE implementation
# The functions are adapted from Edom Moges

nbins <- 11
lower_bd <- 0.25
upper_bd <- 0.42

# This function is to deal with outliers before discretization of continuous data
# Assuming the upper and lower boundaries are always provided for simplification
# Input include:
# The time series need to be discritized: var
# The total number of bins: nbins
# The lower boundary for folding the first bin: lower_bd
# The upper boundary for folding the last bin: upper_bd
# Output1: the number of points in each bin after accounting for the outliers: counts
# Output2: the edges for the bins: breaks
histogram <- function(var,nbins,lower_bd,upper_bd){
  # Get # of outliers on the two sides
  lower_count <- sum(var <= lower_bd,na.rm=TRUE)
  upper_count <- sum(var >= upper_bd,na.rm=TRUE)
  # Filtered data
  filtered <- var[var > lower_bd & var < upper_bd]
  # Get the breaking points for bins
  breaks <- seq(from = lower_bd,to = upper_bd,length.out = nbins + 1)
  # Discretize the data, left-closed bins
  h <- hist(filtered,breaks = breaks,plot=FALSE,right = FALSE)
  # Add outlier counts to the first and last bins
  h$counts[1] <- h$counts[1] + lower_count
  h$counts[length(h$counts)] <- h$counts[length(h$counts)] + upper_count
  return(list(counts = h$counts,
              breaks = breaks))
}

# This function is to discretize continuous data into bins based on bin edges
# The outliers are put into the first/last bins
# Input include:
# The time series need to be discritized: var
# The edges of the bins: binEdges
# The lower boundary for folding the first bin: lower_bd
# The upper boundary for folding the last bin: upper_bd
# Output: returns the bins each point belong to
digitize <- function(var,binEdges,lower_bd,upper_bd){
  # Get the number of bins, this is necessary when adjusting for zero
  nbins <- length(binEdges) - 1
  # Discretize the data, left-closed
  bin_id <- cut(var,breaks = binEdges,include.lowest = TRUE,right = FALSE,labels = FALSE)
  # Assign outliers to the first and last bins
  bin_id[var <= lower_bd] <- 1
  bin_id[var >= upper_bd] <- nbins
  return(bin_id)
}

# This function calculates the upper and lower boundary based on quantile
# Input include:
# The time series need to be discritized: var
# The lower boundary for folding the first bin: lower_qt as in quantile
# The upper boundary for folding the last bin: upper_qt as in quantile
# Output: returns the upper and lower boundary as in their original scale
find_bounds <- function(var,lower_qt,upper_qt){
  lower_bd <- quantile(var,probs = lower_qt,na.rm = TRUE)
  upper_bd <- quantile(var,probs = upper_qt,na.rm = TRUE)
  return(list(lower_bd,upper_bd))
}



