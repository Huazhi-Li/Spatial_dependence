# Script:  Event_simluation_parallel.R
# Details: This script generates an event set of sythetic water levels which are spatially dependent within an user defined cluster. 
#          It first calcultes the dependence structure between sites by the conditional multivariate extreme value model of Heffernan and Tawn, using 3 day maxima water levels.
#          The model is carried out by the function mex of the texmex package. This function consists of two steps. First, Generalized Pareto distributions (GPD) are
#          fitted to the upper tails of each of the marginal distributions of the data: the GPD parameters are estimated for each column of the data in turn, independently 
#          of all other columns. Then, the conditional multivariate approach of Heffernan and Tawn is used to model the dependence between variables.
#          The second part is the stochastic event generation.
# Author:  Huazhi Li, huazhi.li@vu.nl
# Date:    2021-27-09



## Initialization ----------------------------------------------------------
rm(list = ls())

# Load packages
library(texmex)
library(xlsx)
library(stringr)


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
station_id <- as.integer(args[3])

var_in <- paste('Output/10000yr_simulation/', zone, '/', 'cluster', as.character(cluster_id), '/Event_set/',  sep="")
load(paste(var_in, "variables.Rdata", sep="")) 
csv_out <- var_in
dir.create(csv_out)
pth_csv <- paste(csv_out, "pth/", sep="")
dir.create(pth_csv)


i <- station_id

## Dependence model ----------------------------------------------------------
# This model contains two step: (1) Calculating marginal distributions at each gauge and
# (2) Fitting the conditional models between sites. The 'mexAll' function fits a collection of GPD
# and conditional dependence models, by calling the function 'mex' which is a wrapper for calls to migpd and mexDependence

# Calculate dependence between sites through fitting a collection of GPD and conditional dependence models


mexList[[i]] <- mex(data=data,which=i,mqu=mqu,dqu=dqu)

# ## Event simulation ----------------------------------------------------------
# Initialize the event set to be simulated
mult <- 20 # Multiply factor enabling for generating a large pool to randomly sample events 
pqu <- equ #99th prediction threshold


if (count_consite[i]!=0) {

  # Simulate events for each (conditioning) site


  # Reject samples where the conditioning site is not the largest out of all locations (on the quantile scale)
  n_CondLargest <- 0
  n <- 0
  while (n_CondLargest<count_consite[i]) {
    event_sim <- predict.mex(mexList[[i]], which = i, pqu = pqu, nsim = mult*count_consite[i])
    n_CondLargest <- length(which(event_sim$data$CondLargest))
    n = n + 1
    if (n>=5) {
      mult = as.integer(mult*count_consite[i]/n_CondLargest*1.5)
      n <- 0
    }
  }

  whichmax <- which(event_sim$data$CondLargest)
  event_set <- event_sim$data$simulated[sample(whichmax,count_consite[i]),]
  rownames(event_set) <- c(1:nrow(event_set))

  write.csv(event_set, file = paste(csv_out, as.character(i),  ".csv", sep=""), row.names = FALSE)
  
  write.csv(event_sim$data$pth, file = paste(pth_csv, as.character(i),  ".csv", sep=""), row.names = FALSE)
} else {
  event_sim <- predict.mex(mexList[[i]], which = i, pqu = pqu, nsim = 10)
  write.csv(event_sim$data$pth, file = paste(pth_csv, as.character(i),  ".csv", sep=""), row.names = FALSE)
}










