# Script:  Threshold_selection.R
# Details: Selecting approriate threshold for 3-day maxima extreme sea levels at each site, using
#          (1) the cross-validation threshold selection of Northrop et al. (2017), which can be used in package threshr
#          (2) the Stability Threshold Method (STM) and Mean Residual Life (MRL) plot developed by Coles (2001), which can be used in function tstab.gpd (mev) and mrl (texmex) respectively
# Author:  Huazhi Li, huazhi.li@vu.nl
# Date:    2021-31-08

# control+shift+c: put into comment

## Initialization ----------------------------------------------------------
# Load pacakges
options(warn=-1)
library(threshr)
library(ncdf4)
library(stringr)

# Set working directory
setwd("C:/Users/hli490/Desktop/Spatial_dependence/")
file_in <- 'Input/data_gtsm_era5/cf_esl_3d_maxima/'


## Bayesian Cross-Validation method ------------------------------------------------------------

station_ids <- c(seq(0,18716, by = 1), 39421, 39482, 39630)
u_best <- integer(length(station_ids))
j = 1

for (i in station_ids)
  {
  # Read 3-day maxima extreme sea levels
  name_s <- paste('gtsm_station', str_pad(i, 5, pad = "0"), '.nc', sep="") 
  loc_s <- paste(file_in, name_s, sep="")
  esl_3d <- ncvar_get(nc_data, "waterlevel")
  
  
  # Set a vector of training thresholds
  u <- quantile(esl_3d, probs = seq(0.1, 0.998, by = 0.01))
  
  cv_res <- ithresh(data = esl_3d, u_vec = u, n = 1000)
  u_best[j] <- summary(cv_res)[4]
  j = j + 1
  }
