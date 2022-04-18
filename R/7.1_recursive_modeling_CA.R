rm(list=ls())
library(forecast)
library(reshape2)
library(ggplot2)
library(ggmap)
#library(maps)
#library(mapdata)
library(mgcv)
library(TSA)
library(dplyr)
library(fda)
#library(sp)
library(gstat)
library(fields)
library(MBA)
#library(Vizumap)
library(nlme)
library(autoimage)
library(FRK)
library(pwr)
library(MTS)
library(tseries)
library(RColorBrewer)
library(dglm)
library(quantreg)
library(rgbif)
library(mgcv)
library(mgcViz)
library(lctools)
#library(spgwr)

## get aspect ratio from lon/lat
map_aspect = function(x, y) {
  x.center <- sum(range(x)) / 2
  y.center <- sum(range(y)) / 2
  x.dist <- ggplot2:::dist_central_angle(x.center + c(-0.5, 0.5), rep(y.center, 2))
  y.dist <- ggplot2:::dist_central_angle(rep(x.center, 2), y.center + c(-0.5, 0.5))
  y.dist / x.dist
}

tas_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/tas-rec-rcp85-mon-44i/"
pr_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/pr-rec-rcp85-mon-44i/"
map_path = "/Users/wenyilin/Documents/R/NA-CORDEX/map/"
code_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/"
res_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/"
within_ca = readRDS(file = paste0(map_path,"within_ca.rds"))
load(paste0(map_path,"within_rec_ca.rdata"))
load(paste0(map_path,"ca_elevation_rec.rdata"))
load("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/reports/ca_geodata.rdata")

### list files in CA
#coords_ca_all = list.files(path = data_path,pattern = "^coords_ca")
tas_hist_ca_all = list.files(path = tas_path,pattern = "^tas_ca_tas.hist")
tas_rcp_ca_all = list.files(path = tas_path,pattern = "^tas_ca_tas.rcp85")
### example analysis for CA
tas_hist_ca = readRDS(file = paste0(tas_path,tas_hist_ca_all[1]))
tas_rcp_ca = readRDS(file = paste0(tas_path,tas_rcp_ca_all[1]))

# find non-ocean area
coords_ca$ind = 1
tmp.coords = merge(ca_elevation, coords_ca, by = c("latitude","longitude"),all.x=TRUE)
land.ind = which(tmp.coords$ind==1 | tmp.coords$elevation_geonames.x>0)
ca_elevation = ca_elevation[land.ind,]
tas_hist_ca = tas_hist_ca[land.ind,]
tas_rcp_ca = tas_rcp_ca[land.ind,]
ca_elevation$ele_std = (ca_elevation$elevation_geonames - min(ca_elevation$elevation_geonames))/(max(ca_elevation$elevation_geonames)-min(ca_elevation$elevation_geonames))
ca_elevation$ele_norm = (ca_elevation$elevation_geonames - mean(ca_elevation$elevation_geonames))/1000

## seasonal tas data
time_dat = data.frame(year = rep(1950:2005,each=12),
                      month = rep(1:12,56),
                      ts.point = 1:(12*56))
gls.dat = data.frame(lat = ca_elevation[,1], lon = ca_elevation[,2],
                     elevation = ca_elevation[,5],
                     tas = tas_hist_ca,
                     loc = as.factor(1:length(ca_elevation[,1])))
gls.dat.long = reshape2::melt(gls.dat,id.vars = c("lon","lat","elevation","loc"))
names(gls.dat.long) = c("lon","lat","elevation","loc","ts.point","tas")
gls.dat.long$ts.point = as.numeric(gsub("\\D", "", gls.dat.long$ts.point))
gls.dat.long = merge(gls.dat.long,time_dat, by.x = "ts.point",all.x = TRUE)
gls.dat.summary = as.data.frame(gls.dat.long %>% 
                                  filter(month %in% c(1,2,3)) %>%
                                  group_by(lon,lat,elevation,loc,year) %>%
                                  dplyr::summarise(winter = mean(tas))
)
gls.dat.summary$spring = as.data.frame(gls.dat.long %>% 
                                         filter(month %in% c(4,5,6)) %>%
                                         group_by(lon,lat,elevation,loc,year) %>%
                                         dplyr::summarise(spring = mean(tas))
)$spring

gls.dat.summary$summer = as.data.frame(gls.dat.long %>% 
                                         filter(month %in% c(7,8,9)) %>%
                                         group_by(lon,lat,elevation,loc,year) %>%
                                         dplyr::summarise(summer = mean(tas))
)$summer

gls.dat.summary$fall = as.data.frame(gls.dat.long %>% 
                                       filter(month %in% c(10,11,12)) %>%
                                       group_by(lon,lat,elevation,loc,year) %>%
                                       dplyr::summarise(fall = mean(tas))
)$fall
gls.tas.summary = gls.dat.summary[order(gls.dat.summary$loc,gls.dat.summary$year),]

## future data
time_rcp = data.frame(year = rep(2006:2100,each=12),
                      month = rep(1:12,95),
                      ts.point = 1:(12*95))  ## future projection
gls.dat = data.frame(lat = ca_elevation[,1], lon = ca_elevation[,2],
                     elevation = ca_elevation$ele_norm,
                     tas = tas_rcp_ca,
                     loc = as.factor(1:length(ca_elevation[,1])))
gls.dat.long = reshape2::melt(gls.dat,id.vars = c("lon","lat","elevation","loc"))
names(gls.dat.long) = c("lon","lat","elevation","loc","ts.point","tas")
gls.dat.long$ts.point = as.numeric(gsub("\\D", "", gls.dat.long$ts.point))
gls.dat.long = merge(gls.dat.long,time_rcp, by.x = "ts.point",all.x = TRUE)
gls.dat.summary = as.data.frame(gls.dat.long %>% 
                                  filter(month %in% c(1,2,3)) %>%
                                  group_by(lon,lat,elevation,loc,year) %>%
                                  dplyr::summarise(winter = mean(tas))
)
gls.dat.summary$spring = as.data.frame(gls.dat.long %>% 
                                         filter(month %in% c(4,5,6)) %>%
                                         group_by(lon,lat,elevation,loc,year) %>%
                                         dplyr::summarise(spring = mean(tas))
)$spring

gls.dat.summary$summer = as.data.frame(gls.dat.long %>% 
                                         filter(month %in% c(7,8,9)) %>%
                                         group_by(lon,lat,elevation,loc,year) %>%
                                         dplyr::summarise(summer = mean(tas))
)$summer

gls.dat.summary$fall = as.data.frame(gls.dat.long %>% 
                                       filter(month %in% c(10,11,12)) %>%
                                       group_by(lon,lat,elevation,loc,year) %>%
                                       dplyr::summarise(fall = mean(tas))
)$fall
gls.tas.rcp = gls.dat.summary[order(gls.dat.summary$loc,gls.dat.summary$year),]

## color map
cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
neg_cmap = colorRampPalette(rev(brewer.pal(11, "RdBu"))[1:6]) ## all blue
pos_cmap = colorRampPalette(rev((brewer.pal(11, "RdBu"))[1:6])) ## all red

## exploratory plots
season.var = c("winter","spring","summer","fall" )
loc.names = c("SF","SD","Yosemite","Death Valley")
select.loc = data.frame(sf = c(37.75,-122.25), sd = c(32.75,-117.25),
                        yo = c(37.75,-119.25), dv = c(36.25,-116.75))
select.ind = c(100, 1, 106, 78)
select.loc.col = c("red","green","blue","darkgrey")
season.var = c("winter","spring","summer","fall")
sub.loc = list(x = c(-122.25, -117.25, -119.25, -116.75), 
               y = c(37.75,32.75,37.75,36.25),
               labels = loc.names)

n.season = length(season.var)
n.loc = length(unique(gls.tas.summary$loc))
time.yr = (range(gls.tas.summary$year)[1] : range(gls.tas.summary$year)[2])
time.rcp.yr = (range(gls.tas.rcp$year)[1] : range(gls.tas.rcp$year)[2])
#time.pts = c(time.yr - mean(time.yr))/sd(time.yr)
time.pts = c(time.yr - mean(time.yr))/10
time.rcp.pts = c(time.rcp.yr - mean(time.rcp.yr))/10
n.time = length(time.pts)
n.rcp.time = length(time.rcp.pts)

gls.tas.summary$time.pts = time.pts
gls.tas.rcp$time.pts = time.rcp.pts
gls.tas.all = rbind(gls.tas.summary,gls.tas.rcp)

gls.pr.summary$time.pts = time.pts
gls.pr.rcp$time.pts = time.rcp.pts
gls.pr.all = rbind(gls.pr.summary,gls.pr.rcp)

## spatial seasonal trend
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          gls.tas.summary[gls.tas.summary$year==1992,season.var],
          interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Average seasonal temperature (C) in 1992",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          points = sub.loc,
          points.args = list(pch = 20, col = "white"),
          text = sub.loc,
          text.args = list(pos = 3, col = "darkgreen",cex=1.1),
          zlim = c(-10,32),
          col=cmap(200),
          size = c(1, 4), lratio = 0.35,
          legend = "vertical")

## temporal seasonal trend
par(mfrow = c(2,2), oma=c(1, 1, 2,5),  mar = c(2, 2, 2, 2))
for(s in 1:length(select.ind))
{
  tas.dat = ts(gls.tas.summary[gls.tas.summary$loc==select.ind[s],season.var[1]],start = 1950)
  plot(tas.dat, ylim = range(gls.tas.summary[,season.var]), main = loc.names[s],
       xlab = "Year", ylab = "Average temperature")
  for(i in 2: length(season.var))
  {
    tas.dat = ts(gls.tas.summary[gls.tas.summary$loc==select.ind[s],season.var[i]],start = 1950)
    lines(tas.dat, ylim = range(gls.tas.summary[,season.var]),col=i)
  }
}
legend(par('usr')[2], par('usr')[4]*2, bty='n', xpd=NA,cex = 1.1,
       legend=season.var,lty = 1, col=1:length(season.var))
mtext("Average seasonal temperature (C) at four typical locations", outer = TRUE, cex = 1.2)

### combined with future
par(mfrow = c(2,2), oma=c(1, 1, 2,5),  mar = c(2, 2, 2, 2))
for(s in 1:length(select.ind))
{
  tas.dat = ts(gls.tas.summary[gls.tas.summary$loc==select.ind[s],season.var[1]],start = 1950)
  plot(tas.dat, xlim = c(1950,2100),ylim = range(gls.tas.all[,season.var]), main = loc.names[s],
       xlab = "Year", ylab = "Average temperature")
  tas.dat = ts(gls.tas.rcp[gls.tas.rcp$loc==select.ind[s],season.var[1]],start = 2006)
  lines(tas.dat, lty=3)
  for(i in 2: length(season.var))
  {
    tas.dat = ts(gls.tas.summary[gls.tas.summary$loc==select.ind[s],season.var[i]],start = 1950)
    lines(tas.dat, ylim = range(gls.tas.summary[,season.var]),col=i)
    tas.dat = ts(gls.tas.rcp[gls.tas.rcp$loc==select.ind[s],season.var[i]],start = 2006)
    lines(tas.dat, lty=3,col=i)
  }
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=season.var,lty = 1, col=1:length(season.var))
mtext("Average seasonal temperature (C) at three typical locations", outer = TRUE, cex = 1.2)

## map of elevation
par(oma = c(1, 1, 1, 1), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          ca_elevation$elevation_geonames,
          interp.args = list(no.X = 200, no.Y = 200),
          common.legend = FALSE,
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = c("elevation (m) in CA"),
          points = sub.loc,
          points.args = list(pch = 20, col = "white"),
          text = sub.loc,
          text.args = list(pos = 3, col = "white",cex=1.1),
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          #size = c(2,3), lratio = 0.4,
          #col.args = list(cmap1(100),cmap2(100))``
          #col=cmap1(100),
          #zlim = list(c(10,35),c(0,3500)),
          legend = "vertical")

#############################################
a####### Slope in non-weighted models ############
#############################################

## Historical
source(paste0(code_path,"varx_fixed.R"))
fixed_beta1 = rbind(matrix(1,nrow = 2,ncol = 4),
                    c(0,0,0,0),
                    c(0,0,0,0),
                    c(0,0,0,0),
                    c(0,0,0,0))
fixed_beta2 = rbind(matrix(1,nrow = 2,ncol = 4),
                    c(1,0,0,0),
                    c(0,1,0,0),
                    c(0,0,1,0),
                    c(0,0,0,1))

beta_est = beta_se = matrix(0, nrow = n.loc, ncol = n.season)
b0_est = matrix(0, nrow = n.loc, ncol = n.season)
season_coef = season_se = matrix(0, nrow = n.loc, ncol = n.season)
sresids = resids = array(NA, dim = c(n.loc, n.time-1, n.season))
aics = rep(0, n.loc)
aics.lm = rep(0, n.loc)
coef.cov = array(NA, dim = c(n.loc, 8 , 8))
resid_cov = array(NA, dim = c(n.loc, n.season, n.season))
for(i in 1:n.loc)
{
  tas_exp = gls.tas.summary[gls.tas.summary$loc==i,season.var]
  x_multi = cbind(time.pts[-1],tas_exp$fall[1:(n.time-1)],
                  tas_exp$winter[2:n.time],tas_exp$spring[2:n.time],
                  tas_exp$summer[2:n.time])
  x_multi2 = cbind(time.pts[-1],
                  tas_exp$winter[1:(n.time-1)],tas_exp$spring[1:(n.time-1)],
                  tas_exp$summer[1:(n.time-1)],tas_exp$fall[1:(n.time-1)])
  m_ar0_cor2 = varx_fixed(tas_exp[-1,],p=0,xt = x_multi,fixed = fixed_beta2)
  m_ar0_cor1 = varx_fixed(tas_exp[-1,],p=0,xt = x_multi2,fixed = fixed_beta2)
  beta_est[i,] = m_ar0_cor2$beta[,1]
  b0_est[i,] = m_ar0_cor2$coef[1,]
  beta_se[i,] = m_ar0_cor2$se.beta[,1]
  season_coef[i,] = c(m_ar0_cor2$beta[1,2],m_ar0_cor2$beta[2,3],
                      m_ar0_cor2$beta[3,4],
                      m_ar0_cor2$beta[4,5])
  season_se[i,] = c(m_ar0_cor2$se.beta[1,2],m_ar0_cor2$se.beta[2,3],
                    m_ar0_cor2$se.beta[3,4],
                    m_ar0_cor2$se.beta[4,5])
  resids[i,,] = m_ar0_cor2$residuals
  aics[i] = m_ar0_cor2$aic
  aics.lm[i] = m_ar0_cor1$aic
  resid_cov[i,,] = m_ar0_cor2$Sigma
  #coef.cov[i,,] = m_ar0_cor2$coef.cov[[1]]
  #resids_c = center_apply(m_ar0_cor2$residuals)
  sresids[i,,] = t((chol(solve(m_ar0_cor2$Sigma)))%*% t(m_ar0_cor2$residuals))
}

beta_ik = array(0, dim = c(n.loc, n.season, n.season))
b_ik = array(0, dim = c(n.loc, n.season, n.season))
a_i0k = array(0, dim = c(n.loc, n.season, n.season))
a_i1 = array(0, dim = c(n.loc, n.season))
sigma_taui = array(0, dim = c(n.loc,n.season))
se2_taui = array(0, dim = c(n.loc,n.season))
k0 = c(1,2,3,4)
k1 = c(4,1,2,3)
k2 = c(3,4,1,2)
k3 = c(2,3,4,1)
k_ind = cbind(k0,k1,k2,k3)
delta_t = 0.1
grad.beta = grad.b = array(0, dim = c(n.loc,4))

for(i in 1:n.loc)
{
  beta_ik[i,1,] = beta_est[i,k0]
  beta_ik[i,2,] = beta_est[i,k1]
  beta_ik[i,3,] = beta_est[i,k2]
  beta_ik[i,4,] = beta_est[i,k3]
  b_ik[i,1,] = season_coef[i,k0]
  b_ik[i,2,] = season_coef[i,k1]
  b_ik[i,3,] = season_coef[i,k2]
  b_ik[i,4,] = season_coef[i,k3]
  a_i0k[i,,] = rbind(rep(1,4), b_ik[i,1,],
                     b_ik[i,1,] * b_ik[i,2,],
                     b_ik[i,1,] * b_ik[i,2,] * b_ik[i,3,])
  
  ## standard deviation calculator
  cov_tau_tmp = rep(0, n.season)
  #resids_diff = resids[i,2:(n.time-1),] - resids[i,1:(n.time-2),]
  for(ind_1 in 1:n.season)
  {
    for(ind_2 in 1:n.season)
    {
      cov_tau_tmp = cov_tau_tmp + 
        a_i0k[i,ind_1,]*a_i0k[i,ind_2,] * diag(crossprod(resids[i,,k_ind[,ind_1]],resids[i,,k_ind[,ind_2]])/(n.time-1))
      #print(cov_tau_tmp)
    }
  }
  sigma_taui[i,] = cov_tau_tmp
  
  ## standard error calculator
  grad.beta1 = a_i0k[i,k_ind[1,],1]
  grad.beta2 = a_i0k[i,k_ind[2,],2]
  grad.beta3 = a_i0k[i,k_ind[3,],3]
  grad.beta4 = a_i0k[i,k_ind[4,],4]
  
  b.t = season_coef[i,]
  beta.t = beta_est[i,]
  grad.b1 = c(beta.t[4] + beta.t[3] * b.t[4] + beta.t[2] * b.t[4] * b.t[3],
              0, beta.t[2] * b.t[1] * b.t[4],
              beta.t[3] * b.t[1] + beta.t[2] * b.t[1] * b.t[3])
  grad.b2 = c(beta.t[4] * b.t[2] + beta.t[3] * b.t[2] * b.t[4],
              beta.t[1] + beta.t[4] * b.t[1] + beta.t[3] * b.t[1] * b.t[4],
              0, beta.t[3] * b.t[2] * b.t[1])
  grad.b3 = c(beta.t[4] * b.t[3] * b.t[2],
              beta.t[1] * b.t[3] + beta.t[4] * b.t[3] * b.t[1],
              beta.t[2] + beta.t[1] * b.t[2] + beta.t[4] * b.t[2] * b.t[1],0)
  grad.b4 = c(0, beta.t[1] * b.t[4] * b.t[3],
              beta.t[2] * b.t[4] + beta.t[1] * b.t[4] * b.t[2],
              beta.t[3] + beta.t[2] * b.t[3] + beta.t[1] * b.t[3] * b.t[2])
  
  grad.mat = cbind(c(grad.beta1, grad.b1),c(grad.beta2, grad.b2),
                   c(grad.beta3, grad.b3),c(grad.beta4, grad.b4))
  cov.t = diag(c(beta_se[i,]^2,season_se[i,]^2))
  
  se2_taui[i,] = diag(t(grad.mat) %*% cov.t %*% grad.mat)
}

a_i1 = apply(b_ik, c(1,3), prod)
beta_tau = apply(beta_ik * a_i0k, c(1,3), sum) / (1 - a_i1)
sigma_tau = sqrt(se2_taui /(1-a_i1^2))

beta_lm = beta_se_lm = array(0, dim = c(n.loc, 4))
for(i in 1:n.loc)
{
  tas_exp = gls.tas.summary[gls.tas.summary$loc==i,season.var]
  for(j in 1:n.season)
  {
    fit.t = lm(tas_exp[,j] ~ time.pts)
    beta_lm[i,j] = coef(fit.t)[2]
    beta_se_lm[i,j] = summary(fit.t)$coefficients[2,2]
  }
}

## Fit with future data
beta_rcp_est = beta_rcp_se = matrix(0, nrow = n.loc, ncol = n.season)
b0_rcp_est = matrix(0, nrow = n.loc, ncol = n.season)
season_rcp_coef = season_rcp_se = matrix(0, nrow = n.loc, ncol = n.season)
sresids_rcp = resids_rcp = array(NA, dim = c(n.loc, n.rcp.time-1, n.season))
aics_rcp = rep(0, n.loc)
for(i in 1:n.loc)
{
  tas_exp = gls.tas.rcp[gls.tas.rcp$loc==i,season.var]
  x_multi = cbind(time.rcp.pts[-1],tas_exp$fall[1:(n.rcp.time-1)],
                  tas_exp$winter[2:n.rcp.time],tas_exp$spring[2:n.rcp.time],
                  tas_exp$summer[2:n.rcp.time])
  m_ar0_cor2 = varx_fixed(tas_exp[-1,],p=0,xt = x_multi,fixed = fixed_beta2)
  beta_rcp_est[i,] = m_ar0_cor2$beta[,1]
  b0_rcp_est[i,] = m_ar0_cor2$coef[1,]
  beta_rcp_se[i,] = m_ar0_cor2$se.beta[,1]
  season_rcp_coef[i,] = c(m_ar0_cor2$beta[1,2],m_ar0_cor2$beta[2,3],
                          m_ar0_cor2$beta[3,4],
                          m_ar0_cor2$beta[4,5])
  season_rcp_se[i,] = c(m_ar0_cor2$se.beta[1,2],m_ar0_cor2$se.beta[2,3],
                        m_ar0_cor2$se.beta[3,4],
                        m_ar0_cor2$se.beta[4,5])
  resids_rcp[i,,] = m_ar0_cor2$residuals
  aics_rcp[i] = m_ar0_cor2$aic
  #resids_c = center_apply(m_ar0_cor2$residuals)
  sresids_rcp[i,,] = t((chol(solve(m_ar0_cor2$Sigma)))%*% t(m_ar0_cor2$residuals))
}

beta_rcp_ik = array(0, dim = c(n.loc, 4, 4))
b_rcp_ik = array(0, dim = c(n.loc, 4, 4))
a_rcp_i0k = array(0, dim = c(n.loc, 4, 4))
a_rcp_i1 = array(0, dim = c(n.loc, 4))
sigma_rcp_taui = array(0, dim = c(n.loc,n.season))
se2_rcp_taui = array(0, dim = c(n.loc,n.season))
for(i in 1:n.loc)
{
  beta_rcp_ik[i,1,] = beta_rcp_est[i,]
  beta_rcp_ik[i,2,] = beta_rcp_est[i,c(4,1,2,3)]
  beta_rcp_ik[i,3,] = beta_rcp_est[i,c(3,4,1,2)]
  beta_rcp_ik[i,4,] = beta_rcp_est[i,c(2,3,4,1)]
  b_rcp_ik[i,1,] = season_rcp_coef[i,]
  b_rcp_ik[i,2,] = season_rcp_coef[i,c(4,1,2,3)]
  b_rcp_ik[i,3,] = season_rcp_coef[i,c(3,4,1,2)]
  b_rcp_ik[i,4,] = season_rcp_coef[i,c(2,3,4,1)]
  a_rcp_i0k[i,,] = rbind(rep(1,4), b_rcp_ik[i,1,],
                         b_rcp_ik[i,1,] * b_rcp_ik[i,2,],
                         b_rcp_ik[i,1,] * b_rcp_ik[i,2,] * b_rcp_ik[i,3,])
  
  ## standard deviation calculator
  cov_tau_tmp = rep(0, n.season)
  #resids_diff = resids[i,2:(n.time-1),] - resids[i,1:(n.time-2),]
  for(ind_1 in 1:n.season)
  {
    for(ind_2 in 1:n.season)
    {
      cov_tau_tmp = cov_tau_tmp + 
        a_rcp_i0k[i,ind_1,]*a_rcp_i0k[i,ind_2,] * diag(crossprod(resids_rcp[i,,k_ind[,ind_1]],resids_rcp[i,,k_ind[,ind_2]])/(n.time-1))
      #print(cov_tau_tmp)
    }
  }
  sigma_rcp_taui[i,] = cov_tau_tmp
  
  ## standard error calculator
  grad.beta1 = a_rcp_i0k[i,k_ind[1,],1]
  grad.beta2 = a_rcp_i0k[i,k_ind[2,],2]
  grad.beta3 = a_rcp_i0k[i,k_ind[3,],3]
  grad.beta4 = a_rcp_i0k[i,k_ind[4,],4]
  
  b.t = season_rcp_coef[i,]
  beta.t = beta_rcp_est[i,]
  grad.b1 = c(beta.t[4] + beta.t[3] * b.t[4] + beta.t[2] * b.t[4] * b.t[3],
              0, beta.t[2] * b.t[1] * b.t[4],
              beta.t[3] * b.t[1] + beta.t[2] * b.t[1] * b.t[3])
  grad.b2 = c(beta.t[4] * b.t[2] + beta.t[3] * b.t[2] * b.t[4],
              beta.t[1] + beta.t[4] * b.t[1] + beta.t[3] * b.t[1] * b.t[4],
              0, beta.t[3] * b.t[2] * b.t[1])
  grad.b3 = c(beta.t[4] * b.t[3] * b.t[2],
              beta.t[1] * b.t[3] + beta.t[4] * b.t[3] * b.t[1],
              beta.t[2] + beta.t[1] * b.t[2] + beta.t[4] * b.t[2] * b.t[1],0)
  grad.b4 = c(0, beta.t[1] * b.t[4] * b.t[3],
              beta.t[2] * b.t[4] + beta.t[1] * b.t[4] * b.t[2],
              beta.t[3] + beta.t[2] * b.t[3] + beta.t[1] * b.t[3] * b.t[2])
  
  grad.mat = cbind(c(grad.beta1, grad.b1),c(grad.beta2, grad.b2),
                   c(grad.beta3, grad.b3),c(grad.beta4, grad.b4))
  cov.t = diag(c(beta_rcp_se[i,]^2,season_rcp_se[i,]^2))
  
  se2_rcp_taui[i,] = diag(t(grad.mat) %*% cov.t %*% grad.mat)
}

a_rcp_i1 = apply(b_rcp_ik, c(1,3), prod)
beta_rcp_tau = apply(beta_rcp_ik * a_rcp_i0k, c(1,3), sum) / (1 - a_rcp_i1)
sigma_rcp_tau = sqrt(se2_rcp_taui /(1-a_rcp_i1^2))

beta_rcp_lm =beta_rcp_se_lm = array(0, dim = c(n.loc, 4))
for(i in 1:n.loc)
{
  tas_exp = gls.tas.rcp[gls.tas.rcp$loc==i,season.var]
  for(j in 1:n.season)
  {
    #ind_test = (nrow(tas_exp)-20):nrow(tas_exp)
    fit.t = lm(tas_exp[,j] ~ time.rcp.pts)
    beta_rcp_lm[i,j] = coef(fit.t)[2]
    beta_rcp_se_lm[i,j] = summary(fit.t)$coefficients[2,2]
  }
}

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          beta_rcp_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated time trends (C/decade) of four seasons (adjusted with season effects)",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=pos_cmap(200),
          zlim = c(0,1),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          beta_rcp_tau,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons (adjusted with season effects)",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=pos_cmap(200),
          zlim = c(0,1),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          beta_w_rcp_tau,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons (adjusted with season effects)",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=pos_cmap(200),
          zlim = c(0,1),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          beta_rcp_lm,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons (y~t)",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=pos_cmap(200),
          zlim = c(0,1),
          legend = "vertical")

## test convergence of the slope
lags = seq(1,length(time.rcp.pts)-1, by=1)
beta_rcp_lm = array(0, dim = c(n.loc, 4))
for(i in 1:n.loc)
{
  tas_exp = gls.tas.rcp[gls.tas.rcp$loc==i,season.var]
  for(j in 1:n.season)
  {
    for(lag_t in lags)
    {
      
    }
    ind_test = (nrow(tas_exp)-lag_t):nrow(tas_exp)
    fit.t = lm(tas_exp[ind_test,j] ~ time.rcp.pts[ind_test])
    beta_rcp_lm[i,j] = coef(fit.t)[2]
  }
}

#############################################
####### Slope in weighted models ############
#############################################

## Historical
load("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/summary_results_20211101/ca_tas_coef.rdata")
beta_w_ik = array(0, dim = c(n.loc, n.season, n.season))
b_w_ik = array(0, dim = c(n.loc, n.season, n.season))
a_w_i0k = array(0, dim = c(n.loc, n.season, n.season))
a_w_i1 = array(0, dim = c(n.loc, n.season))
sigma_w_taui = array(0, dim = c(n.loc,n.season))
se2_w_taui = array(0, dim = c(n.loc,n.season))
k0 = c(1,2,3,4)
k1 = c(4,1,2,3)
k2 = c(3,4,1,2)
k3 = c(2,3,4,1)
k_ind = cbind(k0,k1,k2,k3)
delta_t = 0.1
grad.beta.w = grad.b.w = array(0, dim = c(n.loc,4))

w_beta_est = t(ca_tas_coef$ca_tas_w$coef[2,,])
w_beta_se = t(ca_tas_coef$ca_tas_w$se.coef[2,,])
w_sea_coef = cbind(c(ca_tas_coef$ca_tas_w$coef[4,1,]),
                   c(ca_tas_coef$ca_tas_w$coef[5,2,]),
                   c(ca_tas_coef$ca_tas_w$coef[6,3,]),
                   c(ca_tas_coef$ca_tas_w$coef[7,4,]))
w_sea_se = cbind(c(ca_tas_coef$ca_tas_w$se.coef[4,1,]),
                   c(ca_tas_coef$ca_tas_w$se.coef[5,2,]),
                   c(ca_tas_coef$ca_tas_w$se.coef[6,3,]),
                   c(ca_tas_coef$ca_tas_w$se.coef[7,4,]))
w_resids = ca_tas_coef$ca_tas_w$residuals
for(i in 1:n.loc)
{
  beta_w_ik[i,1,] = w_beta_est[i,k0]
  beta_w_ik[i,2,] = w_beta_est[i,k1]
  beta_w_ik[i,3,] = w_beta_est[i,k2]
  beta_w_ik[i,4,] = w_beta_est[i,k3]
  b_w_ik[i,1,] = w_sea_coef[i,k0]
  b_w_ik[i,2,] = w_sea_coef[i,k1]
  b_w_ik[i,3,] = w_sea_coef[i,k2]
  b_w_ik[i,4,] = w_sea_coef[i,k3]
  a_w_i0k[i,,] = rbind(rep(1,4), b_w_ik[i,1,],
                     b_w_ik[i,1,] * b_w_ik[i,2,],
                     b_w_ik[i,1,] * b_w_ik[i,2,] * b_w_ik[i,3,])
  
  ## standard deviation calculator
  cov_tau_tmp = rep(0, n.season)
  #resids_diff = resids[i,2:(n.time-1),] - resids[i,1:(n.time-2),]
  for(ind_1 in 1:n.season)
  {
    for(ind_2 in 1:n.season)
    {
      cov_tau_tmp = cov_tau_tmp + 
        a_w_i0k[i,ind_1,]*a_w_i0k[i,ind_2,] * diag(crossprod(w_resids[,k_ind[,ind_1],i],w_resids[,k_ind[,ind_2],i])/(n.time-1))
      #print(cov_tau_tmp)
    }
  }
  sigma_w_taui[i,] = cov_tau_tmp
  
  ## standard error calculator
  grad.beta1 = a_w_i0k[i,k_ind[1,],1]
  grad.beta2 = a_w_i0k[i,k_ind[2,],2]
  grad.beta3 = a_w_i0k[i,k_ind[3,],3]
  grad.beta4 = a_w_i0k[i,k_ind[4,],4]
  
  b.t = w_sea_coef[i,]
  beta.t = w_beta_est[i,]
  grad.b1 = c(beta.t[4] + beta.t[3] * b.t[4] + beta.t[2] * b.t[4] * b.t[3],
              0, beta.t[2] * b.t[1] * b.t[4],
              beta.t[3] * b.t[1] + beta.t[2] * b.t[1] * b.t[3])
  grad.b2 = c(beta.t[4] * b.t[2] + beta.t[3] * b.t[2] * b.t[4],
              beta.t[1] + beta.t[4] * b.t[1] + beta.t[3] * b.t[1] * b.t[4],
              0, beta.t[3] * b.t[2] * b.t[1])
  grad.b3 = c(beta.t[4] * b.t[3] * b.t[2],
              beta.t[1] * b.t[3] + beta.t[4] * b.t[3] * b.t[1],
              beta.t[2] + beta.t[1] * b.t[2] + beta.t[4] * b.t[2] * b.t[1],0)
  grad.b4 = c(0, beta.t[1] * b.t[4] * b.t[3],
              beta.t[2] * b.t[4] + beta.t[1] * b.t[4] * b.t[2],
              beta.t[3] + beta.t[2] * b.t[3] + beta.t[1] * b.t[3] * b.t[2])
  
  grad.mat = cbind(c(grad.beta1, grad.b1),c(grad.beta2, grad.b2),
                   c(grad.beta3, grad.b3),c(grad.beta4, grad.b4))
  cov.t = diag(c(w_beta_se[i,]^2,w_sea_se[i,]^2))
  
  se2_w_taui[i,] = diag(t(grad.mat) %*% cov.t %*% grad.mat)
}

a_w_i1 = apply(b_w_ik, c(1,3), prod)
beta_w_tau = apply(beta_w_ik * a_w_i0k, c(1,3), sum) / (1 - a_w_i1)
sigma_w_tau = sqrt(se2_w_taui /(1-a_w_i1^2))

## Future
beta_w_rcp_ik = array(0, dim = c(n.loc, n.season, n.season))
b_w_rcp_ik = array(0, dim = c(n.loc, n.season, n.season))
a_w_rcp_i0k = array(0, dim = c(n.loc, n.season, n.season))
a_w_rcp_i1 = array(0, dim = c(n.loc, n.season))
sigma_w_rcp_taui = array(0, dim = c(n.loc,n.season))
se2_w_rcp_taui = array(0, dim = c(n.loc,n.season))
k0 = c(1,2,3,4)
k1 = c(4,1,2,3)
k2 = c(3,4,1,2)
k3 = c(2,3,4,1)
k_ind = cbind(k0,k1,k2,k3)
delta_t = 0.1

w_rcp_beta_est = t(ca_tas_coef$ca_tas_pre$coef[2,,])
w_rcp_beta_se = t(ca_tas_coef$ca_tas_pre$se.coef[2,,])
w_rcp_sea_coef = cbind(c(ca_tas_coef$ca_tas_pre$coef[4,1,]),
                   c(ca_tas_coef$ca_tas_pre$coef[5,2,]),
                   c(ca_tas_coef$ca_tas_pre$coef[6,3,]),
                   c(ca_tas_coef$ca_tas_pre$coef[7,4,]))
w_rcp_sea_se = cbind(c(ca_tas_coef$ca_tas_pre$se.coef[4,1,]),
                 c(ca_tas_coef$ca_tas_pre$se.coef[5,2,]),
                 c(ca_tas_coef$ca_tas_pre$se.coef[6,3,]),
                 c(ca_tas_coef$ca_tas_pre$se.coef[7,4,]))
w_rcp_resids = ca_tas_coef$ca_tas_pre$residuals
for(i in 1:n.loc)
{
  beta_w_rcp_ik[i,1,] = w_rcp_beta_est[i,k0]
  beta_w_rcp_ik[i,2,] = w_rcp_beta_est[i,k1]
  beta_w_rcp_ik[i,3,] = w_rcp_beta_est[i,k2]
  beta_w_rcp_ik[i,4,] = w_rcp_beta_est[i,k3]
  b_w_rcp_ik[i,1,] = w_rcp_sea_coef[i,k0]
  b_w_rcp_ik[i,2,] = w_rcp_sea_coef[i,k1]
  b_w_rcp_ik[i,3,] = w_rcp_sea_coef[i,k2]
  b_w_rcp_ik[i,4,] = w_rcp_sea_coef[i,k3]
  a_w_rcp_i0k[i,,] = rbind(rep(1,4), b_w_rcp_ik[i,1,],
                       b_w_rcp_ik[i,1,] * b_w_rcp_ik[i,2,],
                       b_w_rcp_ik[i,1,] * b_w_rcp_ik[i,2,] * b_w_rcp_ik[i,3,])
  
  ## standard deviation calculator
  cov_tau_tmp = rep(0, n.season)
  #resids_diff = resids[i,2:(n.time-1),] - resids[i,1:(n.time-2),]
  for(ind_1 in 1:n.season)
  {
    for(ind_2 in 1:n.season)
    {
      cov_tau_tmp = cov_tau_tmp + 
        a_w_rcp_i0k[i,ind_1,]*a_w_rcp_i0k[i,ind_2,] * diag(crossprod(w_rcp_resids[,k_ind[,ind_1],i],w_rcp_resids[,k_ind[,ind_2],i])/(n.time-1))
      #print(cov_tau_tmp)
    }
  }
  sigma_w_rcp_taui[i,] = cov_tau_tmp
  
  ## standard error calculator
  grad.beta1 = a_w_rcp_i0k[i,k_ind[1,],1]
  grad.beta2 = a_w_rcp_i0k[i,k_ind[2,],2]
  grad.beta3 = a_w_rcp_i0k[i,k_ind[3,],3]
  grad.beta4 = a_w_rcp_i0k[i,k_ind[4,],4]
  
  b.t = w_rcp_sea_coef[i,]
  beta.t = w_rcp_beta_est[i,]
  grad.b1 = c(beta.t[4] + beta.t[3] * b.t[4] + beta.t[2] * b.t[4] * b.t[3],
              0, beta.t[2] * b.t[1] * b.t[4],
              beta.t[3] * b.t[1] + beta.t[2] * b.t[1] * b.t[3])
  grad.b2 = c(beta.t[4] * b.t[2] + beta.t[3] * b.t[2] * b.t[4],
              beta.t[1] + beta.t[4] * b.t[1] + beta.t[3] * b.t[1] * b.t[4],
              0, beta.t[3] * b.t[2] * b.t[1])
  grad.b3 = c(beta.t[4] * b.t[3] * b.t[2],
              beta.t[1] * b.t[3] + beta.t[4] * b.t[3] * b.t[1],
              beta.t[2] + beta.t[1] * b.t[2] + beta.t[4] * b.t[2] * b.t[1],0)
  grad.b4 = c(0, beta.t[1] * b.t[4] * b.t[3],
              beta.t[2] * b.t[4] + beta.t[1] * b.t[4] * b.t[2],
              beta.t[3] + beta.t[2] * b.t[3] + beta.t[1] * b.t[3] * b.t[2])
  
  grad.mat = cbind(c(grad.beta1, grad.b1),c(grad.beta2, grad.b2),
                   c(grad.beta3, grad.b3),c(grad.beta4, grad.b4))
  cov.t = diag(c(w_rcp_beta_se[i,]^2,w_rcp_sea_se[i,]^2))
  
  se2_w_rcp_taui[i,] = diag(t(grad.mat) %*% cov.t %*% grad.mat)
}

a_w_rcp_i1 = apply(b_w_rcp_ik, c(1,3), prod)
beta_w_rcp_tau = apply(beta_w_rcp_ik * a_w_rcp_i0k, c(1,3), sum) / (1 - a_w_rcp_i1)
sigma_w_rcp_tau = sqrt(se2_w_rcp_taui /(1-a_w_rcp_i1^2))

#############################################
############# Save Results ##################
#############################################

slope_hist = list(beta_tau = beta_tau, sigma_tau = sigma_tau,
                  beta_w_tau = beta_w_tau, sigma_w_tau = sigma_w_tau,
                  beta_lm = beta_lm, beta_se_lm = beta_se_lm)

slope_rcp = list(beta_rcp_tau = beta_rcp_tau, sigma_rcp_tau = sigma_rcp_tau,
                 beta_w_rcp_tau = beta_w_rcp_tau, sigma_w_rcp_tau = sigma_w_rcp_tau,
                 beta_rcp_lm = beta_rcp_lm, beta_rcp_se_lm = beta_rcp_se_lm)

slope_new = list(slope_hist = slope_hist, slope_rcp = slope_rcp)
save(slope_new, file =  "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/summary_results_20211101/ca_slope_v_new.rdata")

#############################################
############ Plots & comparison #############
#############################################
load("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/summary_results_20211101/ca_slope_v_new.rdata")

## Historical
slope_est = slope_new$slope_hist$beta_tau
slope_se = slope_new$slope_hist$sigma_tau
z_slope_score = slope_est/slope_se

w_slope_est = slope_new$slope_hist$beta_w_tau
w_slope_se = slope_new$slope_hist$sigma_w_tau
zw_slope_score = w_slope_est/w_slope_se

lm_slope_est = slope_new$slope_hist$beta_lm
lm_slope_se = slope_new$slope_hist$beta_se_lm
zlm_slope_score = lm_slope_est/lm_slope_se

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          lm_slope_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of slopes in simple regression of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,1),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          slope_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes in multivariate regression (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,1),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          w_slope_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes in weighted multivariate regression (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,1),
          legend = "vertical")

## RCP
slope_est = slope_new$slope_rcp$beta_rcp_tau
slope_se = slope_new$slope_rcp$sigma_rcp_tau
z_slope_score = slope_est/slope_se

w_slope_est = slope_new$slope_rcp$beta_w_rcp_tau
w_slope_se = slope_new$slope_rcp$sigma_w_rcp_tau
zw_slope_score = w_slope_est/w_slope_se

lm_slope_est = slope_new$slope_rcp$beta_rcp_lm
lm_slope_se = slope_new$slope_rcp[[6]]
zlm_slope_score = lm_slope_est/lm_slope_se

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          lm_slope_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes in simple regression of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0.4,1),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          slope_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes in multivariate regression (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0.4,1),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          w_slope_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes in weighted multivariate regression (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0.4,1),
          legend = "vertical")


