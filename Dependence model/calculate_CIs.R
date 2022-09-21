# Script:  calculate_CIs.R
# Details: This script calculates the confidence intervals (5th, 95th) for return levels for each station. A range of return periods are considered: 1, 2, 5, 25, 50, 
#          100, 250, 500, 1000.
# Author:  Huazhi Li, huazhi.li@vu.nl
# Date:    2022-27-08


## Initialization ----------------------------------------------------------
rm(list = ls())

# Load packages
library(xlsx)
library(stringr)
library(ncdf4)
library(lmomco)

# Set working directory
setwd('/gpfs/work4/0/FWC2/Spatial_dependence/')
wl_in <- 'Input/data_gtsm_era5/cf_esl_3d_maxima/'
thr_in <- 'Dependence model/Input/thr_all.csv'
CI_out <- 'Dependence model/Confidence_intervals/'

## Calculate CIs ----------------------------------------------------------
thr_all <- read.csv(file =  thr_in, check.names = FALSE) #read thresholds
rp <- c(1,2,5,10,25,50,100,250,500,1000) # return years
t <- 40 # reanalysis data length, years


# set up CI variables
lower_CI <-  matrix(NA, nrow = 1, ncol = length(rp)) # lower CI, 5th percentile
colnames(lower_CI) <- rp # set column index to rps
return_level <- lower_CI # true return level
upper_CI <- lower_CI # upper CI, 95th percentile

# get the station id
args <- commandArgs(trailingOnly=TRUE)
i <- as.integer(args[1]) # the number of stations to be calculated
station_id <- thr_all$station[i]

# read water level time series
name_s <- paste('gtsm_station', str_pad(station_id, 5, pad = "0"), sep="") 
loc_s <- paste(wl_in, name_s, '.nc', sep="")
nc_data <- nc_open(loc_s)
wl_3d <- ncvar_get(nc_data, "waterlevel")
nc_close(nc_data)

# get extremes
thr <- thr_all$mqu[i] # threshold
n <- length(which(wl_3d>thr)) # number of extremes
lambda <- n/t
q <- 1-1/(rp*lambda) # non-exceedance probabilities for rps 
extremes <- wl_3d[which(wl_3d>thr)]

# fit to GPD
lmoments <-  lmoms(extremes)
gpa_para <- pargpa(lmoments)

# calculate CIs
CI <- genci.simple(f=q, gpa_para, n=100, nsim=1000)
lower_CI[1,] <- CI$lwr 
return_level[1,] <- CI$true
upper_CI[1,] <- CI$upr

## Save CIs ----------------------------------------------------------
write.csv(lower_CI, file = paste(CI_out, 'lower_CI/', name_s, ".csv", sep=""), row.names = FALSE)
write.csv(return_level, file = paste(CI_out, 'return_level/', name_s, ".csv", sep=""), row.names = FALSE)
write.csv(upper_CI, file = paste(CI_out, 'upper_CI/', name_s, ".csv", sep=""), row.names = FALSE)



