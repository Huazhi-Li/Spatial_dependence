# Script:  Event_simluation.R
# Details: To be filled
# Author:  Huazhi Li, huazhi.li@vu.nl
# Date:    2021-27-09


## test: NWA 12 clusters. the 5th cluster Canada Quubec

## Initialization ----------------------------------------------------------
library(texmex)
library(xlsx)
library(stringr)
library(ncdf4)

# Set up
setwd('C:/Users/hli490/Desktop/Spatial_dependence/H&T model/Dependence model/')
source('mexTransform.R')
source('predict.mex.R')
source('u2gpd.R')
source('revTransform.R')
data_in <- 'nwa1.esl_3d.Rdata'
load(data_in)
data <- nwa1.esl_3d


## Calculate the annual extreme event counts ----------------------------------------------------------
equ <- .99    #threshold over which the event is considered to be extreme
yr_obs = 40   #the length of observational data (yrs)

# Create data sets of annual event counts for each station and the whole region
n_event <- matrix(0, nrow=yr_obs, ncol=ncol(data)) # Number of extreme events at each station
event_total <- matrix(0, nrow=yr_obs, ncol=1) # Annual event counts for the region
colnames(n_event) <- colnames(data)
t <- 1 

for (i in 1:yr_obs) {
  event_cal <- rep(0, nrow(data))
  if ((i%%4)==0 ) {
    n = 121
  } else {
    n = 122
  }
  for (j in 1:ncol(data)) {
    n_event[i,j] = length(which(data[t:(t+n-1),j]>quantile(data[,j],equ)))
    event_cal[which(data[t:(t+n-1),j]>quantile(data[,j],equ))] = 1
  }
  t = t+n
  event_total[i,1] = sum(event_cal)
}

## Dependence model ----------------------------------------------------------
# This model contains two step: (1) Calculating marginal distributions at each gauge and 
# (2) Fitting the conditional models between sites. The 'mexAll' function fits a collection of GPD 
# and conditional dependence models, by calling the function 'mex' which is a wrapper for calls to migpd and mexDependence 

mqu = .98 # marginal distribution quantile threshold (determined by the cross-validation method)
dqu = .99 # dependence quantile threshold

# Calculate dependence between sites through fitting a collection of GPD and conditional dependence models
mexList <- mexAll(data, mqu = mqu, dqu = rep(dqu, ncol(data)))

 
## Event simulation ----------------------------------------------------------
yr_sim = 1000 # Years to be simulated
n_obs = matrix(0, nrow=1, ncol=ncol(data)) # The number of observed events given ith column is the most extreme (conditioning site)
n_sim = matrix(0, nrow=1, ncol=ncol(data)) # The number of events to be simulated given ith column is the most extreme (conditioning site)
colnames(n_sim) <- colnames(data)
colnames(n_obs) <- colnames(data)

# Calculate the relative likelihood
data.extremes <- data
data.extremes[,] = 0 # Initialized with 0 (means not extreme)

for (i in 1:ncol(data)) {
  data.extremes[which(data[,i]>quantile(data[,i],equ)),i] <- 1   
  dataLaplace <- mexTransform(mexList[[i]]$margins, margins = mexList[[i]]$dependence$margins, method = "mixture")$transformed
}

Laplace.extremes <- dataLaplace*data.extremes
index.extremes <- which(rowSums(Laplace.extremes)>0)

for (i in 1:ncol(data)){
  n_obs[1,i] <- length(which(apply(Laplace.extremes[index.extremes,],1,which.max)==i)) 
  n_sim[1,i] <- n_obs[1,i]/yr_obs*yr_sim
}

# Simulations
event_sim <- vector('list',ncol(data))
n<-rep(NA,ncol(data))

for (i in 1:ncol(data)){
  event_sim[[i]] <- predict.mex(mexList[[i]], which = i, pqu = .99, nsim = n_sim[1,i])
}

event_transformed <- event_sim[[1]]$data$transformed
event_simulated <- event_sim[[1]]$data$simulated
for(i in 2:ncol(data)){
  event_transformed<-rbind(event_transformed,event_sim[[i]]$data$transformed)
  event_simulated<-rbind(event_simulated,event_sim[[i]]$data$simulated)
}
