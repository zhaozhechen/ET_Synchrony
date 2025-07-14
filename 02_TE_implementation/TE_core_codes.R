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
# Output: returns the bins each point belongs to
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

# This function adjusts zero values in the data when discretizing continuous data, true zero values are put in the first bin
# While the non-zero values are put into the rest n-1 bins, and shift 1 bin to the right
# Input are:
# TS data to be processed (var)
# number of total bins for descretization: nbins
# a tiny value for edge extension: ths (default 1e-6)
# lower and upper bd for the outliers
# Output is: returns the bins each point belongs to, after adjusting for zero values
ZeroAdjustment <- function(var,nbins,ths = 1e-6,lower_bd,upper_bd){
  # Initialize a vector to store output bin index
  bin_idx <- rep(NA,length(var))
  # Get index for zero and non-zero values
  zero_idx <- which(var == 0)
  nonzero_idx <- which(var !=0 & !is.na(var))
  # Get non-zero values
  nonzero_values <- var[nonzero_idx]
  # Get histogram info for nonzero values, accounting for outliers
  h <- histogram(nonzero_values,nbins = nbins - 1,lower_bd,upper_bd)
  bin_edges <- h$breaks
  # Expand both edges a little
  bin_edges[1] <- bin_edges[1] - ths
  bin_edges[length(bin_edges)] <- bin_edges[length(bin_edges)] + ths
  # Discretize nonzero values into bins
  bin_id_nonzero <- digitize(nonzero_values,bin_edges,lower_bd,upper_bd)
  # Shift all these nonzero values 1 bin to the right
  bin_idx[nonzero_idx] <- bin_id_nonzero + 1
  # All 0 are put in the first bin
  bin_idx[zero_idx] <- 1
  return(bin_idx)
}

# This function calculates joint bin counts for 3D matrix
# Input includes:
# a 3-D matrix: M (Xlagged,Yt,Yt-1)
# number of total bins for descretization: nbins
# lower and upper bd for the outliers (here they are vectors for the three columns)
# ZFlag: A logical vector of length 3, indicating which column needs zero-adjustment
# Output1: A vector of joint bin counts: N (length nbins^3)
# Output2: corr: correlation between the first two columns (Xlagged and Yt)
joint_entropy3D <- function(M,nbins,lower_bd,upper_bd,ZFlag){
  # Only keep rows that have no NA
  M <- M[complete.cases(M),]
  # Note: use 10e-4 here to be consistent with the python version
  ths <- 10e-4
  # Discretize each column with or without zero-adjustment
  bin_list <- vector("list",3)
  # Loop over each of the columns
  for(i in 1:3){
    # If zero-adjustment is needed
    if(ZFlag[i]){
      # Get bin index
      bin_list[[i]] <- ZeroAdjustment(M[,i],nbins,ths = ths,lower_bd[i],upper_bd[i])
    }else{
      # If no zero-adjustment is needed
      h <- histogram(M[,i],nbins,lower_bd[i],upper_bd[i])
      binEdges <- h$breaks
      # Expand the two edges
      binEdges[1] <- binEdges[1] - ths
      binEdges[length(binEdges)] <- binEdges[length(binEdges)] + ths
      # Get bin index
      bin_list[[i]] <- digitize(M[,i],binEdges,lower_bd[i],upper_bd[i])
    }
  }
  # Convert 3D bin_idx to 1D bin_idx
  joint_bin <- (bin_list[[1]]-1)*nbins^2 + (bin_list[[2]]-1)*nbins + bin_list[[3]]
  # Count number of obs in each bin (joint counts)
  N <- tabulate(joint_bin,nbins^3)
  # Also calculates correlation between the first and second column (Xlagged and Yt)
  corr <- cor(M[,1],M[,2])
  return(list(N=N,corr = corr))
}

# This function calculates Shannon entropy from counts
# This count could be from any dimensions, but needs to be a vector
# Input: count
# Output: Shannon entropy (Unit:bits)
cal_entropy <- function(counts){
  # Only keep positive counts
  counts <- counts[counts>0]
  # Convert to probs
  probs <- counts/sum(counts)
  # Calculate Shannon entropy
  H <- -sum(probs * log2(probs))
  return(H)
}

# This function calculates 1D entropy of X and Y, mutual information (MI(X,Y)), and TE(X->Y)
# Input include:
# a vector of joint bin counts, which is the result from joint_entropy3D
# number of total bins for descretization: nbins
# Output includes:
# Hx, Hy, MI(X,Y), and TE(X->Y)
# MI(X,Y) = H(X) + H(Y) - H(X,Y)
# TE(x->Y) = H(Yt,Yt-1) + H(Yt-1,Xt-lag) - H(Yt-1) - H(Yt,Yt-1,Xt-lag)
cal_info_metrics_3D <- function(N,nbins){
  # If N is not full length, add 0 to the end
  if(length(N) < nbins^3){
    N <- c(N,rep(0,nbins^3 - length(N)))
  }
  # Convert N to 3D array: Xlagged, Yt, and Yt-1
  N3 <- array(N,dim = c(nbins,nbins,nbins))
  
  # 1D marginal counts
  Mxt <- apply(N3, 1, sum)  # sum over dims 2 and 3
  Myt <- apply(N3, 2, sum)  # sum over dims 1 and 3
  Myt_1 <- apply(N3, 3, sum)  # sum over dims 1 and 2
  
  # 2D marginal counts
  Mxtyt <- apply(N3, c(1,2), sum)  # sum over dim 3
  Mytyt_1 <- apply(N3, c(2,3), sum)  # sum over dim 1
  Myt_1xt <- apply(N3, c(1,3), sum)  # sum over dim 2
  
  # Calculate entropies
  Hxt <- cal_entropy(Mxt)
  Hyt <- cal_entropy(Myt)
  Hyt_1 <- cal_entropy(Myt_1)
  Hxtyt <- cal_entropy(Mxtyt)
  Hytyt_1 <- cal_entropy(Mytyt_1)
  Hyt_1xt <- cal_entropy(Myt_1xt)
  Hxtytyt_1 <- cal_entropy(N)
  
  # Calculate MI and TE
  MI <- Hxt + Hyt - Hxtyt
  TE <- Hytyt_1 + Hyt_1xt - Hyt_1 - Hxtytyt_1
  return(list(Hxt = Hxt,Hyt = Hyt, MI = MI, TE = TE))
}

# This function is to shift TS input and generate a shifted matrix
# Input include:
# The original TS of the source: X
# The original TS of the sink: Y
# Lag for X: lag
# The output is a matrix of:[Xlag,Yt,Yt-1]
Lag_Data <- function(X,Y,lag){
  n <- length(X)
  # Shift the TS
  x_lag <- X[1:(n-lag)]
  yt <- Y[(lag+1):n]
  yt_1 <- Y[lag:(n-1)]
  # Combine them
  M <- cbind(x_lag,yt,yt_1)
  return(M)
}

# This function is to shuffle matrix, for calculation of critical TE values
# It shuffles each column in the matrix, while preserving the locations of NA
# Input is the matrix to shuffle: M
# Output: the shuffled matrix
shuffle_matrix <- function(M){
  # Initialize a blank matrix
  M0 <- matrix(NA,nrow=nrow(M),ncol=ncol(M))
  # Loop over the columns
  for(i in 1:ncol(M)){
    col_values <- M[,i]
    non_na_idx <- which(!is.na(col_values))
    # Only shuffle non-NA values
    shuffled_values <- sample(col_values[non_na_idx])
    # Put them pack
    M0[non_na_idx,i] <- shuffled_values
  }
  return(M0)
}

# This function shuffles each column of the input matrix
# Then apply joint_entropy_3D to the shuffled matrix
# Input includes:
# a 3-D matrix: M (Xlagged,Yt,Yt-1)
# number of total bins for descretization: nbins
# lower and upper bd for the outliers (here they are vectors for the three columns)
# ZFlag: A logical vector of length 3, indicating which column needs zero-adjustment
# Output1: A vector of joint bin counts: N (length nbins^3)
# Output2: corr: correlation between the first two columns (Xlagged and Yt)
joint_3D_shuffle <- function(M,nbins,lower_bd,upper_bd,ZFlag){
  # Shuffle the matrix
  Ms <- shuffle_matrix(M)
  # Get joint counts
  results <- joint_entropy3D(Ms,nbins,lower_bd,upper_bd,ZFlag)
  return(results)
}

# This function is to calculate critical 





