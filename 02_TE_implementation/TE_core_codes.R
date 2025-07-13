# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# This code includes core functions for TE implementation
# The functions are adapted from Edom Moges

# This function is to deal with outliers before discretization of continuous data
# Input include:
# The time series need to be discritized: var
# The lower boundary for folding the first bin: lower_qt (as a quantile)
# The upper boundary for folding the last bin: upper_qt (as a quantile)
# Output: TS after dealing with the outliers

