rm(list=ls())
library(reshape2)
library(ncdf4)
library(sf)

data_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/raw/pre_tas/"
setwd(data_path)
all_files = list.files()
nco = ncdf4::nc_open(all_files[1])
lon = ncdf4::ncvar_get(nco, varid = "lon")
lat = ncdf4::ncvar_get(nco, varid = "lat")
ncdf4::nc_close(nco)

# determine which coords are in each state polygon
map_path = "/Users/wenyilin/Documents/R/NA-CORDEX/map/"
within_ca = readRDS(file = paste0(map_path,"within_ca.rds"))
within_co = readRDS(file = paste0(map_path,"within_co.rds"))
within_ks = readRDS(file = paste0(map_path,"within_ks.rds"))
ca_poly = readRDS(file = paste0(map_path,"ca_poly.rds"))
co_poly = readRDS(file = paste0(map_path,"co_poly.rds"))
ks_poly = readRDS(file = paste0(map_path,"ks_poly.rds"))

# create sf coords of lon/lat
dtf = expand.grid(lon, lat)
names(dtf) = c("lon", "lat")
coords = sf::st_as_sf(dtf, coords = c("lon", "lat"),
                      crs = sf::st_crs(ks_poly))
# save data in .rds
save_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/pre_tas/"
for (i in seq_along(all_files)) {
  # open connection and load relevant data
  nco = ncdf4::nc_open(all_files[i])
  message(paste(i, "/", length(all_files), ": ", all_files[i], sep = ""))
  lon = ncdf4::ncvar_get(nco, varid = "lon")
  lat = ncdf4::ncvar_get(nco, varid = "lat")
  coords = expand.grid(lon, lat)
  tas = ncdf4::ncvar_get(nco, varid = "pr")
  # close connection
  ncdf4::nc_close(nco)
  
  # extract relevant data for each state
  #tas = matrix(c(tas) - 273.15, ncol = dim(tas)[3]) # col: time
  tas = matrix(c(tas), ncol = dim(tas)[3]) 
  coords_ca = coords[within_ca,]
  tas_ca = tas[within_ca, ]
  coords_co = coords[within_co,]
  tas_co = tas[within_co, ]
  coords_ks = coords[within_ks,]
  tas_ks = tas[within_ks, ]
  
  # save data
  
  #saveRDS(coords_ca, file = paste0(save_path,"coords_ca_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  saveRDS(tas_ca, file = paste0(save_path,"pr_ca_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  #saveRDS(coords_co, file = paste0("coords_co_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  saveRDS(tas_co, file = paste0(save_path,"pr_co_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  #saveRDS(coords_ks, file = paste0("coords_ks_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  saveRDS(tas_ks, file = paste0(save_path,"pr_ks_", all_files[i], ".rds", sep = ""), compress = "bzip2")
}

### analysis for Kansas
coords_ks = readRDS(file = paste0(data_path,"nacordex_44i_state_data_files/coords_ks_tasmax.hist.CanESM2.CanRCM4.mon.NAM-44i.raw.nc.rds"))
tasmax_ks = readRDS(file = paste0(data_path,"nacordex_44i_state_data_files/tasmax_ks_tasmax.hist.CanESM2.CanRCM4.mon.NAM-44i.raw.nc.rds"))
# some approaches for plotting
sp_tasmax_ks = sp::SpatialPointsDataFrame(coords = coords_ks,
                                          data = as.data.frame(tasmax_ks))
sp::spplot(sp_tasmax_ks, c("V5"))

autoimage::autoimage(x = coords_ks[,1], y = coords_ks[,2], z = tasmax_ks[,1],map = "county",
                     ylab = "Latitude",xlab = "Longitude")

#### Cope set for Kansas
#Confidence level for CoPE sets.
alpha_CoPE = 0.1
#Nominal expected false area ratio for FARE sets.
alpha_FARE = 0.05
#Number of realizations to produce in the MC simulations.
N=1000
lon = unique(coords_ks[,1])
lat = unique(coords_ks[,2])
ta = 1:12
ta = ta - mean(ta)
tb = 1:12
tb = tb - mean(tb)

na = length(ta)
nb = length(tb)
n = na + nb

#Define Design matrix X.
X1 = c(rep(0,na),rep(1,nb))
X2 = rep(1,n)
X3 = c(ta,rep(0,nb))
X4 = c(rep(0,na),tb)
X = cbind(X1,X2,X3,X4)

dat.long = data.frame(lon = coords_ks[,1], lat = coords_ks[,2], val = tasmax_ks[,1])
dat.wide = acast(dat.long, lon~lat, value.var='val')
mask = (!is.na(dat.wide))*1

Y = array(0,c(length(lon),length(lat),n))
for(j in 1:na)
{
  dat.long = data.frame(lon = coords_ks[,1], lat = coords_ks[,2], val = tasmax_ks[,j])
  dat.wide = acast(dat.long, lon~lat, value.var='val')
  Y[,,j] = dat.wide
}

for(j in (1+na):n)
{
  dat.long = data.frame(lon = coords_ks[,1], lat = coords_ks[,2], val = tasmax_ks[,j])
  dat.wide = acast(dat.long, lon~lat, value.var='val')
  Y[,,j] = dat.wide
}

#Apply the algorithm described in the text.
####
#Compute the estimator beta_hat.
M = ginv(t(X) %*% X) %*% t(X)
beta_hat = array(0,c(length(lon),length(lat),4))
for(i in 1:length(lon))
  for(j in 1:length(lat))
    beta_hat[i,j,] = M %*% Y[i,j,]

#Compute the residuals.
R = array(0,c(length(lon),length(lat),n))
for(i in 1:length(lon))
  for(j in 1:length(lat))
    R[i,j,] = Y[i,j,] - X %*% beta_hat[i,j,]

#Compute empirical variance and normalize.
sigma_hat = apply(R,1:2,sd)
R_tilde = R / rep(sigma_hat,n)

#Compute threshold for CoPE sets.
#Ignore values not on the land.
A = (mask==1)
A[is.na(A)] = FALSE
for(i in 1:n) R_tilde[,,i] = R_tilde[,,i] * A

#Apply MC simulations.
P_MC = MC_gauss(R_tilde,N)

#Compute quantile.
a_CoPE=0
while(P_MC(a_CoPE)>alpha_CoPE) a_CoPE = a_CoPE+0.05

#Ignore values not on land.
beta_hat[,,1][mask==0] = -10
#Plot the CoPE and FARE sets.
####
#Setting level.
c=0.2
#Normalized function.
norm_diff = sqrt(n)/2 * (beta_hat[,,1] - c) / sigma_hat 

setEPS()
#Plot CoPE sets.
postscript("CoPE_ks_exp.eps")
par(mar=c(3,3,2,2)+0.1)
image.map(lon,lat,beta_hat[,,1],mask=mask,col = tim.colors(64),horizontal=FALSE,ylab='',xlab='')
drawContour(lon,lat,beta_hat[,,1],c=c,col="purple",lty=1)
drawContour(lon,lat,norm_diff,c=a_CoPE,col="darkred")
drawContour(lon,lat,norm_diff,c=-a_CoPE,col="darkgreen")
dev.off()
