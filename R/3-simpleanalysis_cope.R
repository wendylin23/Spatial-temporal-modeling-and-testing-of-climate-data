library(reshape2)
library(animation)
library(autoimage)

data_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/tas-rcp85-mon-44i/"
map_path = "/Users/wenyilin/Documents/R/NA-CORDEX/map/"
within_ks = readRDS(file = paste0(map_path,"within_ks.rds"))
### list files in Kansas
coords_ks_all = list.files(path = data_path,pattern = "^coords_ks")
tas_hist_ks_all = list.files(path = data_path,pattern = "^tas_ks_tas.hist")
tas_rcp85_ks_all = list.files(path = data_path,pattern = "^tas_ks_tas.rcp85")
### example analysis for Kansas
coords_ks = readRDS(file = paste0(data_path,coords_ks_all[1]))
tas_hist_ks = readRDS(file = paste0(data_path,tas_hist_ks_all[1]))
tas_rcp85_ks = readRDS(file = paste0(data_path,tas_rcp85_ks_all[1]))

### some approaches for plotting
sp_tas_ks = sp::SpatialPointsDataFrame(coords = coords_ks,
                                       data = as.data.frame(tas_hist_ks))
sp::spplot(sp_tas_ks, c("V1"))

autoimage::autoimage(x = coords_ks[,1], y = coords_ks[,2], z = tas_hist_ks[,1],map = "county",
                     ylab = "Latitude",xlab = "Longitude")

## create animation of the temperature
oopt = ani.options(interval = 0.1)
tas_map_plot <- function(x,y,z,z_lim = FALSE)
{
  zlim = c(range(z)[1],range(z)[2])
  if(z_lim == TRUE)
  {
    autoimage::autoimage(x = x, y = y, z = z[,1],map = "county",
                         ylab = "Latitude",xlab = "Longitude",zlim = zlim)
    interval = ani.options('interval')
    for(i in 2:dim(tas_hist_ks)[2])
    {
      autoimage::autoimage(x = x, y = y, z = z[,i],map = "county",
                           ylab = "Latitude",xlab = "Longitude",zlim = zlim)
      ani.pause()
    }
  }
  else
  {
    autoimage::autoimage(x = x, y = y, z = z[,1],map = "county",
                         ylab = "Latitude",xlab = "Longitude")
    interval = ani.options('interval')
    for(i in 2:dim(tas_hist_ks)[2])
    {
      autoimage::autoimage(x = x, y = y, z = z[,i],map = "county",
                           ylab = "Latitude",xlab = "Longitude")
      ani.pause()
    }
  }
}

#saveGIF(tas_map_plot(coords_ks[,1],coords_ks[,2], tas_hist_ks),interval=0.2,movie.name = "tashist_ks.gif")
#saveGIF(tas_map_plot(coords_ks[,1],coords_ks[,2], tas_hist_ks,z_lim = TRUE),interval=0.2,movie.name = "tashist_ks_equal_legend.gif")

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



###############################
##### analysis for CA #########
###############################
coords_ca = readRDS(file = paste0(data_path,"nacordex_44i_state_data_files/coords_ca_tasmax.hist.CanESM2.CanRCM4.mon.NAM-44i.raw.nc.rds"))
tasmax_ca = readRDS(file = paste0(data_path,"nacordex_44i_state_data_files/tasmax_ca_tasmax.hist.CanESM2.CanRCM4.mon.NAM-44i.raw.nc.rds"))
sp_tasmax_ca = sp::SpatialPointsDataFrame(coords = coords_ca,
                                          data = as.data.frame(tasmax_ca))
sp::spplot(sp_tasmax_ca, c("V1"))
autoimage::autoimage(x = coords_ca[,1], y = coords_ca[,2], z = tasmax_ca[,1],map = "county",
                     ylab = "Latitude",xlab = "Longitude")
