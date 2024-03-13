# Script:  AnnualEventCount_ConditioningStation.R
# Details: This script calculates (1) the total number of events to be simulated and (2) the corresponding conditioning stations.
#          It first extracts historic extreme events from the input 3-day maxima water level series. The number of extreme events are fitted to a distribution using 
#          the Kernel Density Distribution for estimating the total number of events to be generated.
#          The second part is estimating the number of events during which the station i is the conditioning station. To do so, we first calculate the probability that 
#          each station is the conditioning station (i.e. the station has the largest Laplace value) using the empirical data. Then, we assign the total number of events
#          to each conditioning site based on the probability distribution.
# Author:  Huazhi Li, huazhi.li@vu.nl
# Date:    2021-27-09


## 

## Initialization ----------------------------------------------------------
rm(list = ls())

# Load packages
library(texmex)
library(xlsx)
library(stringr)
library(fitur)
library(Ake)
library(EnvStats)

# Add functions
setwd('/gpfs/work4/0/FWC2/Spatial_dependence/Dependence_model/')
source('mexTransform.R')
source('predict.mex.R')
source('u2gpd.R')
source('revTransform.R')

# Load input 3d maximum water levels
args <- commandArgs(trailingOnly=TRUE)
zone <- as.character(args[1])
cluster_id <- as.integer(args[2])

thr_in <- paste('Input/', zone, '/thr.csv', sep="")
data_in <- paste('Input/', zone, '/', zone, as.character(cluster_id), '.Rdata', sep="")
data_out0 <- paste('Output/10000yr_simulation/', zone, '/', sep="")
if (!dir.exists(data_out0)){dir.create(data_out0)}
data_out <- paste(data_out0, 'cluster', as.character(cluster_id), '/', sep="")
if (!dir.exists(data_out)){dir.create(data_out)}
load(data_in)
data <- esl_3d


## Dependence model ----------------------------------------------------------
# This model contains two step: (1) Calculating marginal distributions at each gauge and 
# (2) Fitting the conditional models between sites. The 'mexAll' function fits a collection of GPD 
# and conditional dependence models, by calling the function 'mex' which is a wrapper for calls to migpd and mexDependence 
df_thr <- read.csv(file=thr_in,check.names = FALSE) 
index_thr <- which(df_thr$cluster==cluster_id)
thr <- data.frame("station" = df_thr$station[index_thr],
                  "mth" = df_thr$mth[index_thr])

equ <- .99 # quantile threshold to select extreme events
mqu <- sapply(1:ncol(data), function(i)thr$mth[which(thr$station==colnames(data)[i])]/100) # marginal distribution quantile threshold (determined by the cross-validation method), above which the data is fitted to a generalized pareto distribution and below which data is empirical
dqu <- .92 # dependence quantile threshold
# Set up the dependence calculation
mexList <- vector('list', length=ncol(data))

mexList[[1]] <- mex(data=data,which=1,mqu=mqu,dqu=dqu)


## Define extremes, annual event counts, and conditioning site pool ----------------------------------------------------------
## Calculate the event counts to be simulated for each station
yr_sim = 10000 # Years to be simulated
yr_obs = 40   #the length of observational data (yrs)

# Identify the extreme events by applying peaks over threshold

data.extremes <- data
data.extremes[,] = 0 # Initialized with 0 (means not extreme)

for (i in 1:ncol(data)) {
  data.extremes[which(data[,i]>=quantile(data[,i],equ)),i] <- 1   
}

# Convert to Laplace scale
dataLaplace <- mexTransform(mexList[[1]]$margins, margins = mexList[[1]]$dependence$margins, method = "mixture")$transformed
Laplace.extremes <- dataLaplace*data.extremes
index.extremes <- which(rowSums(data.extremes)>0)

# Calculate the number of observed events for each year of the observation
yr_index <- as.integer(index.extremes/(365.25/3))+1
n_event_obs <- sapply(1:yr_obs,function(i)length(which(yr_index==i)))

# Fit to a non-parametric distribution using the Kernel Density Distribution which can handle discrete data (Wansouwé et al., 2015)
bw <- hcvd.fun(n_event_obs, NULL, "bino") # estimate the bandwidth
pmf <- kpmfe.fun(n_event_obs,bw$hcv,"discrete","bino") # calculate the probability mass function

# Sample annual event counts for simulation
n_event_sim <- sample(x = pmf$eval.points, yr_sim, replace = T, prob = pmf$est.fn) # sample from this fitted distribution

# Calculate the pool of conditioning site (i.e. per event which station has the highest Larplace value)
site_pool <- sapply(1:sum(n_event_obs), function(i)which.max(Laplace.extremes[index.extremes[i],]))

# calculate the probability of each site being the conditioning site
prob_site <- sapply(1:ncol(data), function(i)length(which(site_pool==i))/length(index.extremes))

# generate a condition sample for simulations following the multinomial distribution with ns trials and probabilities Pr.Y 2 Ei /= Pr.Y 2 E/ for i 2 (calculated from the empirical data)
count_consite <- rmultinom(1, size = sum(n_event_sim), prob = prob_site)

# save variables
var_out <- paste(data_out, "Event_set/", sep="")

if (!dir.exists(var_out)){dir.create(var_out)}

save(list = ls(), file = paste(var_out, "/variables.Rdata", sep=""))
