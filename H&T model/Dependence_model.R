# Script:  Dependence_model.R
# Details: Calculate the dependence structure between sites by the conditional multivariate extreme value model of Heffernan and Tawn, using 3 day maxima water levels.
#          The model is carried out by the function mex of the texmex package. This function consists of two steps. First, Generalized Pareto distributions (GPD) are
#          fitted to the upper tails of each of the marginal distributions of the data: the GPD parameters are estimated for each column of the data in turn, independently 
#          of all other columns. Then, the conditional multivariate approach of Heffernan and Tawn is used to model the dependence between variables.
# Author:  Huazhi Li, huazhi.li@vu.nl
# Date:    2021-14-09


## test: NWA 12 clusters. the 5th cluster Canada Quubec

## Initialization ----------------------------------------------------------
library(texmex)
library(xlsx)
library(stringr)
library(ncdf4)
library(MultiHazard)

# Set working directory
#data_in <- '/lustre5/0/FWC2/Spatial_dependence/Input/data_gtsm_era5/cf_esl_3d_maxima/'
#cluster_in <- '/lustre5/0/FWC2/Spatial_dependence/Dependence_model/NWA12.xlsx'
#data_out <- '/lustre5/0/FWC2/Spatial_dependence/Dependence_model/'
data_in <- 'C:/Users/hli490/Desktop/Spatial_dependence/Input/data_gtsm_era5/cf_esl_3d_maxima/'
cluster_in <- 'C:/Users/hli490/Desktop/Spatial_dependence/H&T model/Dependence model/NWA12.xlsx'
data_out <- 'C:/Users/hli490/Desktop/Spatial_dependence/H&T model/Dependence model/'

# Load data
data_cluster <- read.xlsx(cluster_in, sheetIndex = 1)
#cluster_id = 5
#num_s = length(which(data_cluster$Cluster==cluster_id))
num_s = 10
nwa1.esl_3d <- matrix(nrow = 4870, ncol = num_s) # declustered and dethrended 3d maxima 
#colnames(nwa1.esl_3d) <- data_cluster$Station[which(data_cluster$Cluster==cluster_id)]
colnames(nwa1.esl_3d) <- data_cluster$Station[1:10]
j = 1
#which(data_cluster$Cluster==cluster_id)
for (i in seq(1,10,1)) {
  # Read 3-day maxima extreme sea levels
  name_s <- paste('gtsm_station', str_pad(data_cluster$Station[i], 5, pad = "0"), sep="") 
  loc_s <- paste(data_in, name_s, '.nc', sep="")
  nc_data <- nc_open(loc_s)
  wl_3d <- ncvar_get(nc_data, "waterlevel")
  nwa1.esl_3d[,j] = wl_3d
  j = j + 1
}

save(nwa1.esl_3d, file =  paste(data_out, 'nwa1.esl_3d.Rdata', sep=""))

## Dependence model ----------------------------------------------------------
# Calculate marginal distributions independently at each station
thr = .98 # quantile threshold calculated using the cross-validation method
res_migpd <- migpd(nwa1.esl_3d,
                 mqu = .98,
                 penalty = "gaussian",
                 maxit = 10000,
                 trace = 0,
                 verbose = FALSE,
                 priorParameters = NULL,
                 cov = "observed",
                 family = gpd)




# Calculate pairwise dependence structure
# select the first column as the conditioning site
res_dep <- mexDependence(res_migpd, 
                         which = colnames(nwa1.esl_3d)[1],
                         dqu = .99,
                         margins = "laplace",
                         constrain = TRUE,
                         v = 10
                         )


res_pre <- predict(res_dep)

# Simulation
dq_u <- rep(.99, dim(nwa1.esl_3d)[2])

res_mexall <- mexAll(nwa1.esl_3d, mqu = .98, dqu = dq_u)

event_sim <- mexMonteCarlo(nSample = 10, mexList = res_mexall, mult = 10)
save(event_sim, file =  paste(data_out, 'event_sim.Rdata', sep="") )
