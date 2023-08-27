rm(list=ls())
library(forecast)
library(reshape2)
library(ggplot2)
library(ggmap)
library(mgcv)
library(TSA)
library(dplyr)
library(fda)
library(gstat)
library(fields)
library(MBA)
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

code_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/R/"
data_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Data" 
res_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Results"
setwd(code_path)
source(paste0(code_path,"gwr.R"))
source(paste0(code_path,"varx_fixed.R"))
source(paste0(code_path,"util.R"))

#####################################
##### Load and prepare data #########
#####################################
map_path = paste0(data_path,'/map/')
within_ca = readRDS(file = paste0(map_path,"within_ca.rds"))
load(paste0(map_path,"within_rec_ca.rdata"))
load(paste0(map_path,"ca_elevation_rec.rdata"))
load(paste0(map_path,"ca_geodata.rdata"))
### example analysis for CA
tas_hist_ca = readRDS(file = paste0(data_path,'/pre_tas/tas_ca_temp.hist.CanESM2.CanRCM4.mon.NAM-44i.raw.nc.rds'))
tas_rcp_ca = readRDS(file = paste0(data_path,'/pre_tas/tas_ca_temp.rcp85.CanESM2.CanRCM4.mon.NAM-44i.raw.nc.rds'))

# find non-ocean area
coords_ca$ind = 1
tmp.coords = merge(ca_elevation, coords_ca, by = c("latitude","longitude"),all.x=TRUE)
land.ind = which(tmp.coords$ind==1 | tmp.coords$elevation_geonames.x>0)
ca_elevation = ca_elevation[land.ind,]
tas_hist_ca = tas_hist_ca[land.ind,]
tas_rcp_ca = tas_rcp_ca[land.ind,]
ca_elevation$ele_std = (ca_elevation$elevation_geonames - min(ca_elevation$elevation_geonames))/(max(ca_elevation$elevation_geonames)-min(ca_elevation$elevation_geonames))
ca_elevation$ele_norm = (ca_elevation$elevation_geonames - mean(ca_elevation$elevation_geonames))/1000
#ca_elevation = ca_elevation[-282,]

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
gls.tas.rcp = gls.tas.rcp[complete.cases(gls.tas.rcp),]

## color map
cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
neg_cmap = colorRampPalette(rev(brewer.pal(11, "RdBu"))[1:6]) ## all blue
pos_cmap = colorRampPalette(rev((brewer.pal(11, "RdBu"))[1:6])) ## all red

## exploratory plots
season.var = c("winter","spring","summer","fall" )
loc.names = c("San Francisco","San Diego","Yosemite","Death Valley")
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
time.pts = c(time.yr - mean(time.yr))/10
time.rcp.pts = c(time.rcp.yr - mean(time.rcp.yr))/10
n.time = length(time.pts)
n.rcp.time = length(time.rcp.pts)

gls.tas.summary$time.pts = time.pts
gls.tas.rcp$time.pts = time.rcp.pts
gls.tas.all = rbind(gls.tas.summary,gls.tas.rcp)


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
          legend = "vertical")

########################################
###### Temperature Model Fitting  ######
########################################
## 1.non-weighted regression
beta_est = beta_se = matrix(0, nrow = n.loc, ncol = n.season)
b0_est = matrix(0, nrow = n.loc, ncol = n.season)
season_coef = season_se = matrix(0, nrow = n.loc, ncol = n.season)
sresids = resids = array(NA, dim = c(n.loc, n.time-1, n.season))
aics = rep(0, n.loc)
for(i in 1:n.loc)
{
  tas_exp = gls.tas.summary[gls.tas.summary$loc==i,season.var]
  x_multi = cbind(time.pts[-1],tas_exp$fall[1:(n.time-1)],
                  tas_exp$winter[2:n.time],tas_exp$spring[2:n.time],
                  tas_exp$summer[2:n.time])
  fixed_beta2 = rbind(matrix(1,nrow = 2,ncol = 4),
                      c(1,0,0,0),
                      c(0,1,0,0),
                      c(0,0,1,0),
                      c(0,0,0,1))
  m_ar0_cor2 = varx_fixed(tas_exp[-1,],p=0,xt = x_multi,fixed = fixed_beta2)
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
  sresids[i,,] = t((chol(solve(m_ar0_cor2$Sigma)))%*% t(m_ar0_cor2$residuals))
}
## non-weight parameter list
ca_tas_nw = list(b0_est = b0_est,
                 beta_est = beta_est,
                 beta_se = beta_se,
                 season_coef = season_coef,
                 season_se = season_se,
                 resids = resids)

## Fit with future data with non-weighted model
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
  fixed_beta2 = rbind(matrix(1,nrow = 2,ncol = 4),
                      c(1,0,0,0),
                      c(0,1,0,0),
                      c(0,0,1,0),
                      c(0,0,0,1))
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
ca_tas_rcp_nw = list(b0_est = b0_est,
                 beta_est = beta_est,
                 beta_se = beta_se,
                 season_coef = season_coef,
                 season_se = season_se,
                 resids = resids)

## 2. test spatial correlation
n.lon = length(unique(gls.tas.summary$lon))
n.lat = length(unique(gls.tas.summary$lat))
simgrid = cbind(as.numeric(as.factor(ca_elevation$longitude)),
                as.numeric(as.factor(ca_elevation$latitude)))
distance = as.matrix(dist(simgrid))
r.distance = as.matrix(dist(simgrid)) * 50 ##0.44 degree = 50km

## correlation linear/quadratic test
cor.se = 1/sqrt(dim(resids)[2]-3)
cor.list = cor.z = cor.test = list()
cor.lm.x1 = cor.lm.x2 = matrix(0, nrow=n.season, ncol=n.loc)
cor.p.x1 = cor.p.x2 = matrix(0, nrow=n.season, ncol=n.loc)
lmtest.cor = matrix(0, nrow=n.season, ncol=n.loc)
for(t in 1:n.season)
{
  cor.list[[t]] = cor(t(resids[,,t]))
  cor.z[[t]] = atanh(cor.list[[t]])
  cor.test[[t]] = pnorm(cor.z[[1]],lower.tail = FALSE)
  for(i in 1:n.loc)
  {
    tmp.lm1 = lm(cor.list[[t]][i,] ~ distance[i,])
    cor.lm.x1[t,i] = summary(tmp.lm1)$coefficients[2,1]
    cor.p.x1[t,i] = round(summary(tmp.lm1)$coefficients[2,4],5)
    tmp.lm2 = lm(cor.list[[t]][i,] ~ distance[i,] + I(distance[i,]^2))
    cor.lm.x2[t,i] = summary(tmp.lm2)$coefficients[3,1]
    cor.p.x2[t,i] = round(summary(tmp.lm2)$coefficients[3,4],5)
    lmtest.cor[t,i] = round(anova(tmp.lm2,tmp.lm1)$`Pr(>F)`[2],5)
  }
}

op = par(mfrow = c(2,2),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(t in 1:n.season)
{
  tmp.cor = as.vector(cor.list[[t]])
  vec.dist = as.vector(distance)
  tmp.dat = data.frame(tmp.cor = tmp.cor, vec.dist = vec.dist)
  plot(vec.dist, tmp.cor,xlab="distance in pixel",ylab="correlation", 
       pch=".", cex=0.8, main = season.var[t],ylim=c(0,1), col="grey")
  tmp.fit = lm(I(tmp.cor-1) ~ 0 + vec.dist, data = tmp.dat)
  abline(1, coef(tmp.fit), lty=2)
  
  for(i in 1:length(select.ind))
  {
    tmp.fit1 = lm(I(cor.list[[t]][select.ind[i],]-1) ~ 0 + distance[select.ind[i],])
    abline(1, coef(tmp.fit1), col=select.loc.col[i])
  }
  if(t == n.season)
    legend("bottomleft", legend = c("overall",loc.names),
           col = c("black", select.loc.col), lty=c(2,1,1,1,1),cex=0.8)
  #tmp.quad = lm(tmp.cor ~ vec.dist + I(vec.dist^2), data = tmp.dat)
  #dist = seq(0, max(vec.dist), by=0.1)
  #pred.cor = predict(tmp.quad,newdata = data.frame(vec.dist = dist))
  #lines(dist, pred.cor, col="blue")
  #legend("topright", legend = c("linear","quadratic"),col = c("red","blue"), lty=1)
}
title(xlab = "Distance in grids",
      ylab = "Residual correlation",
      outer = TRUE, line = 3)
par(op)

apply(cor.p.x1,1,function(x) table(x<0.05)) ## linear term significant?
apply(cor.p.x2,1,function(x) table(x<0.05)) ## quadratic term significant?

## test kernel functions
m.type = c("log_lm","log_quad","lm","quad")
cor.rsquare = matrix(0, nrow = length(m.type), ncol=n.season)
par(mfrow=c(1,2))
for(t in 1:n.season)
{
  tmp.cor = as.vector(cor.list[[t]])
  vec.dist = as.vector(distance)
  dist.seq = seq(0, max(vec.dist), by=0.1)
  tmp.dat = data.frame(tmp.cor = tmp.cor, vec.dist = vec.dist,
                       log.cor = log(tmp.cor),dist2 = vec.dist^2)
  fit_loglm = lm(log.cor ~ -1 + vec.dist, data = tmp.dat)
  fit_logquad = lm(log.cor ~ -1 + dist2, data = tmp.dat)
  fit_lm1 = lm(I(tmp.cor - 1) ~ 0 + vec.dist, data = tmp.dat) 
  #fit_lm2 = lm(I(sqrt(1-tmp.cor)) ~ 0 + vec.dist, data = tmp.dat)
  fit_lm21 = lm(I(1-tmp.cor) ~ 0 + dist2, data = tmp.dat)
  cor.rsquare[,t] = c(summary(fit_loglm)$r.squared, summary(fit_logquad)$r.squared,
                      summary(fit_lm1)$r.squared, summary(fit_lm21)$r.squared)
  plot(tmp.dat$vec.dist, tmp.dat$log.cor,pch=".", cex=0.8,
       main = season.var[t], xlab = "distance",ylab = "log(cor)")
  abline(fit_loglm, col="red")
  pred.cor = predict(fit_logquad,newdata = data.frame(dist2 = dist.seq^2))
  lines(dist.seq, pred.cor, col="blue")

  plot(tmp.dat$vec.dist, tmp.dat$tmp.cor,pch=".", cex=0.8,
       main = season.var[t], xlab = "distance",ylab = "correlation")
  pred.cor = predict(fit_lm1,newdata = data.frame(vec.dist = dist.seq))
  lines(dist.seq, pred.cor+1, col="red")
  pred.cor = predict(fit_lm2,newdata = data.frame(vec.dist = dist.seq))
  # lines(dist.seq, (1-pred.cor^2), col="blue")
  
  #tmp.quad = lm(tmp.cor ~ vec.dist + I(vec.dist^2), data = tmp.dat)
  #dist = seq(0, max(vec.dist), by=0.1)
  #pred.cor = predict(tmp.quad,newdata = data.frame(vec.dist = dist))
  #lines(dist, pred.cor, col="blue")
  #legend("topright", legend = c("linear","quadratic"),col = c("red","blue"), lty=1)
}

rownames(cor.rsquare) = m.type
colnames(cor.rsquare) = season.var
cor.rsquare

## 3 Constructing GWMTSR 
d_bis = seq(1,10,by = 0.5)
w2 = gwr.lm(distance,d_bis[])

### Model selection via AIC
fixed_beta2 = rbind(matrix(1,nrow = 3,ncol = 4),diag(1,4))
fixed_beta3 = rbind(matrix(1,nrow = 3,ncol = 4),diag(0,4))
aic.bisquare.w1 = aic.bisquare.w2 = matrix(0, nrow = 4, ncol = length(d_bis))
for(j in 1:length(d_bis))
{
  print(j)
  w2 = gwr.lm(distance,d_bis[j])
  fit_t2 = varw_fixed(tas_exp,p=0,xt = x_mult,w = w2,fixed = fixed_beta2)
  aic.bisquare.w2[,j] = fit_t2$aic
}

d_gauss = seq(1,10,by = 0.5)
aic.gauss.w3 = aic.gauss.w4 = matrix(0, nrow = 4, ncol = length(d_gauss))
for(j in 1:length(d_gauss))
{
  print(j)
  w3 = gwr.exp(distance, d_gauss[j] )
  w4 = gwr.gauss(distance^2,d_gauss[j])
  fit_t3 = varw_fixed(tas_exp,p=0,xt = x_mult,w = w3,fixed = fixed_beta2)
  fit_t4 = varw_fixed(tas_exp,p=0,xt = x_mult,w = w4,fixed = fixed_beta2)
  aic.gauss.w3[,j] = fit_t3$aic
  aic.gauss.w4[,j] = fit_t4$aic
}## d[6]=0.6

par(mfrow=c(1,1))
plot(d_bis,colSums(aic.bisquare.w1),type="l",xlab = "d",ylab = "AIC")
lines(d_bis, colSums(aic.bisquare.w2),col="red")
lines(d_bis, colSums(aic.gauss.w3),col="blue")
lines(d_bis,colSums(aic.gauss.w4),col="darkgreen")
legend("topright",legend=c("lm","quad","exp","gauss"),
       col = c("red","black","blue","darkgreen"),lty=1)

aic.list = list(aic.bisquare.w1 = aic.bisquare.w1,
                aic.bisquare.w2 = aic.bisquare.w2,
                aic.gauss.w3 = aic.gauss.w3,
                aic.gauss.w4 = aic.gauss.w4)

### Model selection via 10-fold CV
set.seed(2021)
shuffle.ind = sample(dim(tas_exp)[3])
cv.folds = cut(seq(1,dim(tas_exp)[3]),breaks=10,labels=FALSE)
pred.dat = function(a,b)
{
  pred.y = array(0, dim=c(dim(a)[1],dim(b)[2],dim(a)[3]))
  n.test = dim(a)[3]
  for(i in 1:n.test)
  {
    pred.y[,,i] = a[,,i] %*% b[,,i]
  }
  pred.y
}

n.fold = 10
pred.err = matrix(0, nrow =  length(d_bis), ncol = n.fold)
set.seed(2021)
for(i in 1:n.fold)
{
  # for(j in 1:length(d_bis))
  # {
    w2 = gwr.lm(distance,d_bis[j])
    test.ind = shuffle.ind[which(cv.folds==i)]
    train.ind = shuffle.ind[-which(cv.folds==i)]
    train.ind = train.ind[order(train.ind, decreasing = FALSE)]
    if(j==1)
    {
      w_t = w0
      w_t[test.ind,test.ind] = 0
      fit_t2 = varw_fixed(tas_exp,p=0,xt = x_mult[,-3,],w = w_t,fixed = fixed_beta1) ## no elevation
      predict.resi = pred.dat(x_mult[,-3,test.ind], fit_t2$coef[,,test.ind])
    }
    else{
      w_t = w2
      w_t[test.ind,test.ind] = 0
      fit_t2 = varw_fixed(tas_exp,p=0,xt = x_mult,w = w_t,fixed = fixed_beta2)
      predict.resi = pred.dat(x_mult[,,test.ind], fit_t2$coef[,,test.ind])
    }
    pred.err[j,i] = mean(apply(predict.resi - tas_exp[,,test.ind], 3, function(x) sum(x^2)))
  #}## d=1.6
}

## Bisquare kernel d_bis[2]=1.1
fixed_beta2 = rbind(matrix(1,nrow = 3,ncol = 4),diag(1,4))
fixed_beta3 = rbind(matrix(1,nrow = 3,ncol = 4),diag(0,4))
fixed_beta4 = rbind(matrix(1,nrow = 2,ncol = 4),rep(0,4),diag(1,4))
w1 = gwr.bisquare(distance^2,4)
w2 = gwr.lm(distance,2)
w3 = gwr.exp(distance, 1.5)
w4 = gwr.gauss(distance^2,2)

fixed_beta0 = rbind(matrix(1,nrow = 3,ncol = 4),diag(0,4))
fit_m3_2 = varw_fixed(tas_exp,p=0,xt = x_mult,w = w2,fixed = fixed_beta2) ## selected weight matrix
### Fit the model with future data
pred_m3_2 = varw_fixed(tas_rcp,p=0,xt = x_rcp,w = w2,fixed = fixed_beta2)
pred_m3_beta0 = varw_fixed(tas_rcp,p=0,xt = x_rcp,w = w2,fixed = fixed_beta0)
pred_m3_w3_beta0 = varw_fixed(tas_rcp,p=0,xt = x_rcp,w = w2,fixed = fixed_beta0)
pred_m3_w0_beta0 = varw_fixed(tas_rcp,p=0,xt = x_rcp,w = w0,fixed = fixed_beta0)

### Future data Model selection via AIC
d_bis = seq(0.5,5,by = 0.5)
fixed_beta2 = rbind(matrix(1,nrow = 3,ncol = 4),diag(1,4))
fixed_beta3 = rbind(matrix(1,nrow = 3,ncol = 4),diag(0,4))
fixed_beta4 = rbind(matrix(1,nrow = 2,ncol = 4),rep(0,4),diag(1,4))
aic.bisquare.w1 = aic.bisquare.w2 = matrix(0, nrow = 4, ncol = length(d_bis))
for(j in 1:length(d_bis))
{
  w2 = gwr.lm(distance,d_bis[j])
  fit_t2 = varw_fixed(tas_rcp,p=0,xt = x_rcp,w = w2,fixed = fixed_beta2)
  fit_t3 = varw_fixed(tas_rcp,p=0,xt = x_rcp,w = w2,fixed = fixed_beta3)
  fit_t4 = varw_fixed(tas_rcp,p=0,xt = x_rcp,w = w2,fixed = fixed_beta4)
  aic.bisquare.w2[,j] = fit_t2$aic
  print(fit_t2$aic)
}## d=1.6

d_gauss = seq(1,10,by = 0.5)
aic.gauss.w3 = aic.gauss.w4 = matrix(0, nrow = 4, ncol = length(d_gauss))
for(j in 1:length(d_gauss))
{
  print(j)
  w3 = gwr.exp(distance, d_gauss[j] )
  w4 = gwr.gauss(distance^2,d_gauss[j])
  fit_t3 = varw_fixed(tas_exp,p=0,xt = x_mult,w = w3,fixed = fixed_beta2)
  fit_t4 = varw_fixed(tas_exp,p=0,xt = x_mult,w = w4,fixed = fixed_beta2)
  aic.gauss.w3[,j] = fit_t3$aic
  aic.gauss.w4[,j] = fit_t4$aic
}## d[6]=0.6

## Prediction with future data
rcp_pred = pred.dat(x_rcp, pred_m3_2$coef)
rcp_beta0_pred = pred.dat(x_rcp, pred_m3_beta0$coef)
dat_pred = pred.dat(x_mult, fit_m3_2$coef)

## current prediction
par(mfrow = c(2,2), oma=c(1, 1, 2,5),  mar = c(2, 2, 2, 2))
for(s in 1:length(select.ind))
{
  tas.dat = ts(gls.tas.summary[gls.tas.summary$loc==select.ind[s],season.var[1]],start = 1950)
  plot(tas.dat, ylim = range(gls.tas.all[,season.var]), main = loc.names[s],
       xlab = "Year", ylab = "Average temperature")
  tas.dat = ts(dat_pred[,1,select.ind[s]],start = 1950)
  lines(tas.dat, lty = 2)
  for(i in 2: length(season.var))
  {
    tas.dat = ts(gls.tas.summary[gls.tas.summary$loc==select.ind[s],season.var[i]],start = 1950)
    lines(tas.dat,col=i)
    tas.dat = ts(dat_pred[,i,select.ind[s]],start = 1950)
    lines(tas.dat, lty = 2,col=i)
  }
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c(season.var),lty = 1, col=1:length(season.var))
mtext("Recorded seasonal temperature vs. predicted temperature (1950-2005)", outer = TRUE, cex = 1.2)

## future prediction
par(mfrow = c(2,2), oma=c(1, 1, 2,5),  mar = c(2, 2, 2, 2))
for(s in 1:length(select.ind))
{
  tas.dat = ts(gls.tas.rcp[gls.tas.rcp$loc==select.ind[s],season.var[1]],start = 2006)
  plot(tas.dat, ylim = range(gls.tas.all[,season.var]), main = loc.names[s],
       xlab = "Year", ylab = "Average temperature")
  tas.dat = ts(rcp_pred[,1,select.ind[s]],start = 2006)
  lines(tas.dat, lty = 2)
  for(i in 2: length(season.var))
  {
    tas.dat = ts(gls.tas.rcp[gls.tas.rcp$loc==select.ind[s],season.var[i]],start = 2006)
    lines(tas.dat,col=i)
    tas.dat = ts(rcp_pred[,i,select.ind[s]],start = 2006)
    lines(tas.dat, lty = 2,col=i)
  }
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c(season.var),lty = 1, col=1:length(season.var))
mtext("RCM seasonal temperature vs. predicted temperature (2006-2100)", outer = TRUE, cex = 1.2)

