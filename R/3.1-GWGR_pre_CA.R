###########################################################################################
## This code file provides the analysis pipeline for precipitation data in CA
###########################################################################################

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
pr_hist_ca = readRDS(file = paste0(data_path,'/pre_tas/pr_ca_prec.hist.CanESM2.CanRCM4.mon.NAM-44i.raw.nc.rds'))
pr_rcp_ca = readRDS(file = paste0(data_path,'/pre_tas/pr_ca_prec.rcp85.CanESM2.CanRCM4.mon.NAM-44i.raw.nc.rds'))

## find non-ocean area
coords_ca$ind = 1
tmp.coords = merge(ca_elevation, coords_ca, by = c("latitude","longitude"),all.x=TRUE)
land.ind = which(tmp.coords$ind==1 | tmp.coords$elevation_geonames.x>0)
ca_elevation = ca_elevation[land.ind,]
pr_hist_ca = pr_hist_ca[land.ind,]
pr_rcp_ca = pr_rcp_ca[land.ind,]
ca_elevation$ele_std = (ca_elevation$elevation_geonames - min(ca_elevation$elevation_geonames))/(max(ca_elevation$elevation_geonames)-min(ca_elevation$elevation_geonames))
ca_elevation$ele_norm = (ca_elevation$elevation_geonames - mean(ca_elevation$elevation_geonames))/1000

## seasonal tas data
time_dat = data.frame(year = rep(1950:2005,each=12),
                      month = rep(1:12,56),
                      ts.point = 1:(12*56))
gls.dat = data.frame(lat = ca_elevation[,1], lon = ca_elevation[,2],
                     elevation = ca_elevation[,5],
                     tas = pr_hist_ca*30, #*86400*30,
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
gls.pr.summary = gls.dat.summary[order(gls.dat.summary$loc,gls.dat.summary$year),]

## future data
time_rcp = data.frame(year = rep(2006:2100,each=12),
                      month = rep(1:12,95),
                      ts.point = 1:(12*95))  ## future projection
gls.dat = data.frame(lat = ca_elevation[,1], lon = ca_elevation[,2],
                     elevation = ca_elevation$ele_norm,
                     tas = pr_rcp_ca*30, #*86400*30,
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
gls.pr.rcp = gls.dat.summary[order(gls.dat.summary$loc,gls.dat.summary$year),]

## basic definitions
season.var = c("winter","spring","summer","fall" )
loc.names = c("San Francisco","San Diego","Yosemite","Death Valley")
select.loc = data.frame(sf = c(37.75,-122.25), sd = c(32.75,-117.25),
                        yo = c(37.75,-119.25), dv = c(36.25,-116.75))
select.ind = c(126, 1, 132, 87)
select.loc.col = c("red","green","blue","darkgrey")
season.var = c("winter","spring","summer","fall")
sub.loc = list(x = c(-122.25, -117.25, -119.25, -116.75), 
               y = c(37.75,32.75,37.75,36.25),
               labels = loc.names)
n.season = length(season.var)
n.loc = length(unique(gls.pr.summary$loc))
time.yr = (range(gls.pr.summary$year)[1] : range(gls.pr.summary$year)[2])
time.rcp.yr = (range(gls.pr.rcp$year)[1] : range(gls.pr.rcp$year)[2])
#time.pts = c(time.yr - mean(time.yr))/sd(time.yr)
time.pts = c(time.yr - mean(time.yr))/10
time.rcp.pts = c(time.rcp.yr - mean(time.rcp.yr))/10
n.time = length(time.pts)
n.rcp.time = length(time.rcp.pts)

gls.pr.summary$time.pts = time.pts
gls.pr.rcp$time.pts = time.rcp.pts
gls.pr.all = rbind(gls.pr.summary,gls.pr.rcp)
gls.pr.summary[which(is.na(gls.pr.summary[,"winter"])),season.var]=0
gls.pr.all[which(is.na(gls.pr.all[,"winter"])),season.var]=0

## color map
cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
rev.cmap <- colorRampPalette(brewer.pal(11, "RdBu"))
neg_cmap = colorRampPalette(rev(brewer.pal(11, "RdBu"))[1:6]) ## all blue
pos_cmap = colorRampPalette(rev((brewer.pal(11, "RdBu"))[1:6])) ## all red
## color map for data
pr.cmap <- colorRampPalette(brewer.pal(11, "Spectral")[6:11])
tas.cmap <- colorRampPalette(rev(brewer.pal(11, "Spectral"))[1:6])
ele.cmap <- colorRampPalette(terrain.colors(10, alpha = 1))

## exploratory plots
## spatial seasonal trend
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          gls.pr.summary[gls.pr.summary$year==1992,season.var],
          interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Average seasonal precipitation (mm/month) in 1992",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          points = sub.loc,
          points.args = list(pch = 20, col = "white"),
          text = sub.loc,
          text.args = list(pos = 3, col = "darkred",cex=1.1),
          zlim = c(0,180),
          col=pr.cmap(200),
          size = c(2, 2), lratio = 0.2,
          legend = "vertical")

### combined with future
par(mfrow = c(2,2), oma=c(1, 1, 2,1),  mar = c(2, 2, 2, 2))
for(s in 1:length(select.ind))
{
  tas.dat = ts(gls.pr.summary[gls.pr.summary$loc==select.ind[s],season.var[1]],start = 1950)
  plot(tas.dat, xlim = c(1950,2100),ylim = range(gls.pr.all[,season.var]), main = loc.names[s],
       xlab = "Year", ylab = "Average precipitation")
  tas.dat = ts(gls.pr.rcp[gls.pr.rcp$loc==select.ind[s],season.var[1]],start = 2006)
  lines(tas.dat, lty=3)
  for(i in 2: length(season.var))
  {
    tas.dat = ts(gls.pr.summary[gls.pr.summary$loc==select.ind[s],season.var[i]],start = 1950)
    lines(tas.dat,col=i)
    tas.dat = ts(gls.pr.rcp[gls.pr.rcp$loc==select.ind[s],season.var[i]],start = 2006)
    lines(tas.dat, lty=3,col=i)
  }
  if(s==2)
    legend("topright", bty='n', xpd=NA,cex = 0.9,
           legend=season.var,lty = 1, col=1:length(season.var))
}
#legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
#       legend=season.var,lty = 1, col=1:length(season.var))
mtext("Average seasonal precipitation (mm/month) at four typical locations", outer = TRUE, cex = 1.2)

########################################
###### Precipitation Model Fitting  ####
########################################
## 1. non-weighted gamma regression at each location
shape_est = slope_est = slope_se =  p.value = gof.test = matrix(0, nrow = n.loc, ncol = n.season)
resids = array(NA, dim = c(n.loc, n.time, n.season))
sig.level = 0.1
for(i in 1:n.loc)
{
  for(j in 1:n.season)
  {
    tas_exp = gls.pr.summary[gls.pr.summary$loc == i,season.var[j]]
    fit1 = glm(tas_exp ~ time.pts, family=Gamma(link = "log"))
    resids[i,,j] = fit1$residuals
    deviance = fit1$deviance
    p.value[i,j] = pchisq(deviance, df = fit1$df.residual, lower.tail = F)
    gof.test[i,j] = ifelse(p.value[i,j] < sig.level,0,1)
    shape_est[i,j] <- 1/summary(fit1)$dispersion
    scale_est <- 1/fit1$coefficients[1]/shape_est
    slope_est[i,j] <- fit1$coefficients[2]
    slope_se[i,j] <- summary(fit1)$coefficients[2,2]
  }
}

## 2. geographically weighted gamma regression
gwr_gamma = function(b0,b1,b2, logshape) {
  loglik = rep(0, n.loc)
  for(j in 1:n.loc)
  {
    linear_predictor = x[[j]][,1]*b0 + x[[j]][,2]*b1 + x[[j]][,3]*b2
    sum_j = sum((exp(logshape)-1) *log(y[j,]) - log(gamma(exp(logshape))) - 
                  exp(logshape)*(linear_predictor - logshape) - exp(logshape) * y[j,]/exp(linear_predictor))
    loglik[j] = w_i[j] * sum_j
  }
  # rate = shape / mean:
  # sum of negative log likelihoods:
  -sum(loglik)
}

gr_gamma = function(b0,b1, logshape) {
  #n = dim(y)[1]
  #ts = dim(y)[2]
  loglik = rep(0, n.loc)
  for(j in 1:n.loc)
  {
    linear_predictor = x[[j]][,1]*b0 + x[[j]][,2]*b1
    sum_j = sum((exp(logshape)-1) *log(y[j,]) - log(gamma(exp(logshape))) - 
                  exp(logshape)*(linear_predictor - logshape) - exp(logshape) * y[j,]/exp(linear_predictor))
    loglik[j] = w_i[j] * sum_j
  }
  # rate = shape / mean:
  # sum of negative log likelihoods:
  -sum(loglik)
}

gwr_loglik = function(b0,b1,b2, logshape)
{
  loglik = rep(0, n.loc)
  for(j in 1:n.loc)
  {
    linear_predictor = x[[j]][,1]*b0 + x[[j]][,2]*b1 + x[[j]][,3]*b2
    sum_j = sum((exp(logshape)-1) *log(y[j,]) - log(gamma(exp(logshape))) - 
                  exp(logshape)*(linear_predictor - logshape) - exp(logshape) * y[j,]/exp(linear_predictor))
    loglik[j] = sum_j
  }
  sum(loglik)
}

n.lon = length(unique(gls.pr.summary$lon))
n.lat = length(unique(gls.pr.summary$lat))
simgrid = cbind(as.numeric(as.factor(ca_elevation$longitude)),
                as.numeric(as.factor(ca_elevation$latitude)))
distance = as.matrix(dist(simgrid))
r.distance = as.matrix(dist(simgrid)) * 50 ##0.44 degree = 50km

### residual correlation checking
## Correlation between residuals
cor.se = 1/sqrt(dim(resids)[2]-3)
cor.glm = cor.gwr = list()
cor.glm.test = cor.gwr.test = list()
for(t in 1:n.season)
{
  cor.glm[[t]] = cor(t(resids[,,t]))
  cor.glm.test[[t]] = pnorm(atanh(cor.glm[[t]]),lower.tail = FALSE)
}

op = par(mfrow = c(2,2),
         oma = c(5,4,0,0) + 0.1,
         mar = c(0,0,1,1) + 0.1)
for(t in 1:n.season)
{
  tmp.cor = as.vector(cor.glm[[t]])
  vec.dist = as.vector(distance)
  tmp.dat = data.frame(tmp.cor = tmp.cor, vec.dist = vec.dist)
  plot(vec.dist, tmp.cor,xlab="distance in pixel",ylab="correlation", 
       pch=".", cex=0.8, main = season.var[t],ylim=c(0,1), col="grey")
  tmp.fit = lm(I(tmp.cor-1) ~ 0 + vec.dist, data = tmp.dat)
  abline(1, coef(tmp.fit), lty=2)
}
title(xlab = "distance in grids",
      ylab = "residual correlation",
      outer = TRUE, line = 3)
par(op)

## test kernel functions
m.type = c("log_lm","log_quad","lm","quad")
cor.rsquare = matrix(0, nrow = length(m.type), ncol=n.season)
par(mfrow=c(1,2))
for(t in 1:n.season)
{
  tmp.cor = as.vector(cor.glm[[t]])
  vec.dist = as.vector(distance)
  dist.seq = seq(0, max(vec.dist), by=0.1)
  tmp.dat = data.frame(tmp.cor = tmp.cor, vec.dist = vec.dist,
                       log.cor = log(tmp.cor),dist2 = vec.dist^2)
  fit_loglm = lm(log.cor ~ -1 + vec.dist, data = tmp.dat)
  fit_logquad = lm(log.cor ~ -1 + dist2, data = tmp.dat)
  fit_lm1 = lm(I(tmp.cor - 1) ~ 0 + vec.dist, data = tmp.dat) 
  fit_lm2 = lm(I(sqrt(1-tmp.cor)) ~ 0 + vec.dist, data = tmp.dat)
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
  lines(dist.seq, (1-pred.cor^2), col="blue")
  
  #tmp.quad = lm(tmp.cor ~ vec.dist + I(vec.dist^2), data = tmp.dat)
  #dist = seq(0, max(vec.dist), by=0.1)
  #pred.cor = predict(tmp.quad,newdata = data.frame(vec.dist = dist))
  #lines(dist, pred.cor, col="blue")
  #legend("topright", legend = c("linear","quadratic"),col = c("red","blue"), lty=1)
}
rownames(cor.rsquare) = m.type
colnames(cor.rsquare) = season.var
cor.rsquare

## specify weighted matrix options
d_bis = seq(1,10,by = 0.5)
w0 = diag(1, n.loc) # weight matrix with no spatial correlation
w1 = gwr.lm(distance,d_bis[3]) # weight matrix from linear kernel
w3 = gwr.exp(distance, 6) # weight matrix from exponential kernel

## Data (current) prepared for model fitting
x = list()
y.all = list()
for(i in 1:n.loc)
{
  x[[i]] = cbind(rep(1, n.time), time.pts, gls.pr.summary$elevation[gls.pr.summary$loc==i])
}

for(j in 1:n.season)
{
  y.all[[j]] = matrix(NA, nrow = n.loc, ncol = n.time)
  for(i in 1:n.loc)
  {
    y.all[[j]][i,] =  gls.pr.summary[gls.pr.summary$loc == i,season.var[j]]
  }
}

pr_exp = array(0, dim = c(n.time,n.season,n.loc))
pr_rcp = array(0, dim = c(n.rcp.time,n.season,n.loc))
x_mult = array(0, dim = c(n.time,3,n.loc))
x_rcp = array(0, dim = c(n.rcp.time,3,n.loc))
for(i in 1:n.loc)
{
  data.i = gls.pr.summary[gls.pr.summary$loc==i,]
  data.rcp.i = gls.pr.rcp[gls.pr.rcp$loc==i,]
  pr_exp[,,i] = as.matrix(data.i[,season.var])
  pr_rcp[,,i] = as.matrix(data.rcp.i[,season.var])
  x_mult[,,i] =  cbind(rep(1,n.time),
                       time.pts,
                       rep(data.i$elevation[1],n.time))
  x_rcp[,,i] =  cbind(rep(1,n.rcp.time),
                      time.rcp.pts,
                      rep(data.rcp.i$elevation[1],n.rcp.time))
}

## GWGR with weight matrix
gwgr.b0 = gwgr.slope = gwgr.h = matrix(0, nrow = n.loc, ncol=n.season)
gwgr.b0.se = gwgr.slope.se = gwgr.h.se = matrix(0, nrow = n.loc, ncol=n.season)
gwgr.loglik = matrix(0, nrow = n.loc, ncol=n.season)
dat_pred = dat_resid = pr_exp
rcp_pred = pr_rcp
for(i in 1:n.loc)
{
  w_i = w1[i,]
  for(j in 1:n.season)
  {
    y = y.all[[j]]
    m_mle2 = bbmle::mle2(gwr_gamma,optimfun = "BFGS",
                         start = list(b0 = 2,b1 = slope_est[i,j], b2=0.01,
                                      logshape = log(shape_est[i,j])))
    gwgr.b0[i,j] = bbmle::coef(m_mle2)[1]
    gwgr.slope[i,j] = bbmle::coef(m_mle2)[2]
    gwgr.h[i,j] = bbmle::coef(m_mle2)[3]
    gwgr.b0.se[i,j] = sqrt(diag(bbmle::vcov(m_mle2)))[1]
    gwgr.slope.se[i,j] = sqrt(diag(bbmle::vcov(m_mle2)))[2]
    gwgr.h.se[i,j] = sqrt(diag(bbmle::vcov(m_mle2)))[3]
    gwgr.loglik[i,j] = gwr_loglik(bbmle::coef(m_mle2)[1],bbmle::coef(m_mle2)[2],
                                  bbmle::coef(m_mle2)[3],bbmle::coef(m_mle2)[4])
    beta.t = bbmle::coef(m_mle2)[1:3]
    dat_pred[,j,i] = exp(x_mult[,,i]%*%beta.t)[,1]
    dat_resid[,j,i] = y[i,]-dat_pred[,j,i]
    rcp_pred[,j,i] = exp(x_rcp[,,i]%*%beta.t)[,1]
  }
}

res.gwgr = list(gwgr.b0 = gwgr.b0, gwgr.slope = gwgr.slope, gwgr.h=gwgr.h,
                gwgr.b0.se=gwgr.b0.se,gwgr.slope.se=gwgr.slope.se,gwgr.h.se=gwgr.h.se,
                gwgr.loglik = gwgr.loglik,
                dat_pred = dat_pred,dat_resid = dat_resid)

## GWGR without weight matrix (equivalent with GR)
gr.b0 = gr.slope = gr.h = matrix(0, nrow = n.loc, ncol=n.season)
gr.b0.se = gr.slope.se = gr.h.se = matrix(0, nrow = n.loc, ncol=n.season)
gr.loglik = matrix(0, nrow = n.loc, ncol=n.season)
dat_pred = dat_resid = pr_exp
for(i in 1:n.loc)
{
  w_i = w0[i,]
  for(j in 1:n.season)
  {
    y = y.all[[j]]
    m_mle2 = bbmle::mle2(gr_gamma,optimfun = "BFGS",
                         start = list(b0 = 2,b1 = slope_est[i,j],
                                      logshape = log(shape_est[i,j])))
    gr.b0[i,j] = bbmle::coef(m_mle2)[1]
    gr.slope[i,j] = bbmle::coef(m_mle2)[2]
    #gr.h[i,j] = coef(m_mle2)[3]
    gr.b0.se[i,j] = sqrt(diag(bbmle::vcov(m_mle2)))[1]
    gr.slope.se[i,j] = sqrt(diag(bbmle::vcov(m_mle2)))[2]
    #gr.h.se[i,j] = sqrt(diag(bbmle::vcov(m_mle2)))[3]
    gr.loglik[i,j] = gwr_loglik(bbmle::coef(m_mle2)[1],bbmle::coef(m_mle2)[2],
                                0,bbmle::coef(m_mle2)[3])
    beta.t = bbmle::coef(m_mle2)[1:3]
    dat_pred[,j,i] = exp(x_mult[,,i]%*%beta.t)[,1]
    dat_resid[,j,i] = y[i,]-dat_pred[,j,i]
  }
}

res.gr = list(gwgr.b0 = gr.b0, gwgr.slope = gr.slope, gwgr.h=gr.h,
              gwgr.b0.se=gr.b0.se,gwgr.slope.se=gr.slope.se,gwgr.h.se=gr.h.se,
              gwgr.loglik = gr.loglik,
              dat_pred = dat_pred,dat_resid = dat_resid)

### Data (future) prepared for model fitting
x = list()
y.all = list()
for(i in 1:n.loc)
{
  x[[i]] = cbind(rep(1, n.rcp.time), time.rcp.pts, gls.pr.rcp$elevation[gls.pr.rcp$loc==i])
}

for(j in 1:n.season)
{
  y.all[[j]] = matrix(NA, nrow = n.loc, ncol = n.rcp.time)
  for(i in 1:n.loc)
  {
    y.all[[j]][i,] =  gls.pr.rcp[gls.pr.rcp$loc == i,season.var[j]]
  }
}
gwgr.b0 = gwgr.slope = gwgr.h = matrix(0, nrow = n.loc, ncol=n.season)
gwgr.b0.se = gwgr.slope.se = gwgr.h.se = matrix(0, nrow = n.loc, ncol=n.season)
gwgr.loglik = matrix(0, nrow = n.loc, ncol=n.season)
dat_pred = pr_exp
rcp_pred = rcp_resid = pr_rcp
for(i in 1:n.loc)
{
  w_i = w1[i,]
  for(j in 1:n.season)
  {
    y = y.all[[j]]
    m_mle2 = bbmle::mle2(gwr_gamma,optimfun = "BFGS",
                         start = list(b0 = 2,b1 = slope_est[i,j], b2=0.01,
                                      logshape = log(shape_est[i,j])))
    gwgr.b0[i,j] = bbmle::coef(m_mle2)[1]
    gwgr.slope[i,j] = bbmle::coef(m_mle2)[2]
    gwgr.h[i,j] = bbmle::coef(m_mle2)[3]
    gwgr.b0.se[i,j] = sqrt(diag(bbmle::vcov(m_mle2)))[1]
    gwgr.slope.se[i,j] = sqrt(diag(bbmle::vcov(m_mle2)))[2]
    gwgr.h.se[i,j] = sqrt(diag(bbmle::vcov(m_mle2)))[3]
    gwgr.loglik[i,j] = gwr_loglik(bbmle::coef(m_mle2)[1],bbmle::coef(m_mle2)[2],
                                  bbmle::coef(m_mle2)[3],bbmle::coef(m_mle2)[4])
    beta.t = bbmle::coef(m_mle2)[1:3]
    rcp_pred[,j,i] = exp(x_rcp[,,i]%*%beta.t)[,1]
    rcp_resid[,j,i] = y[i,]-rcp_pred[,j,i]
  }
}

future.gwgr = list(gwgr.b0 = gwgr.b0, gwgr.slope = gwgr.slope, gwgr.h=gwgr.h,
                   gwgr.b0.se=gwgr.b0.se,gwgr.slope.se=gwgr.slope.se,gwgr.h.se=gwgr.h.se,
                   gwgr.loglik = gwgr.loglik,
                   dat_pred = rcp_pred,dat_resid = rcp_resid)

## Save data to local directory for figure generation
gwgr.pre.ca = list(res.gr = res.gr,
                   res.gwgr = res.gwgr,
                   future.gwgr = future.gwgr
                   #gwgr.AICc.lm = gwgr.AICc
                   )
#save(gwgr.pre.ca,file = "yourdatapath/ca_pre_coef.rdata")


## current prediction
par(mfrow = c(2,2), oma=c(1, 1, 2,5),  mar = c(2, 2, 2, 2))
for(s in 1:length(select.ind))
{
  tas.dat = ts(gls.pr.summary[gls.pr.summary$loc==select.ind[s],season.var[1]],start = 1950)
  plot(tas.dat, ylim = c(0,300), main = loc.names[s],
       xlab = "Year", ylab = "Average precipitation")
  tas.dat = ts(dat_pred[,1,select.ind[s]],start = 1950)
  lines(tas.dat, lty = 2)
  for(i in 2: length(season.var))
  {
    tas.dat = ts(gls.pr.summary[gls.pr.summary$loc==select.ind[s],season.var[i]],start = 1950)
    lines(tas.dat,col=i)
    tas.dat = ts(dat_pred[,i,select.ind[s]],start = 1950)
    lines(tas.dat, lty = 2,col=i)
  }
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c(season.var),lty = 1, col=1:length(season.var))
mtext("Recorded seasonal precipitation vs. predicted precipitation (1950-2005)", outer = TRUE, cex = 1)

## future prediction
par(mfrow = c(2,2), oma=c(1, 1, 2,5),  mar = c(2, 2, 2, 2))
for(s in 1:length(select.ind))
{
  tas.dat = ts(gls.pr.rcp[gls.pr.rcp$loc==select.ind[s],season.var[1]],start = 2006)
  plot(tas.dat, ylim = range(gls.pr.all[,season.var]), main = loc.names[s],
       xlab = "Year", ylab = "Average precipitation")
  tas.dat = ts(rcp_pred[,1,select.ind[s]],start = 2006)
  lines(tas.dat, lty = 2)
  for(i in 2: length(season.var))
  {
    tas.dat = ts(gls.pr.rcp[gls.pr.rcp$loc==select.ind[s],season.var[i]],start = 2006)
    lines(tas.dat,col=i)
    tas.dat = ts(rcp_pred[,i,select.ind[s]],start = 2006)
    lines(tas.dat, lty = 2,col=i)
  }
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c(season.var),lty = 1, col=1:length(season.var))
mtext("RCM seasonal precipitation vs. predicted precipitation (2006-2100)", outer = TRUE, cex = 1)


## Compute AICc
n.df = 4
gwgr.AICc = apply(gwgr.loglik, 2, function(x) -2*sum(x) + 2*n.df*n.loc + 2*(n.df*n.loc)*(n.df*n.loc+1)/(n.time*n.loc)-n.df*n.loc-1)
n.df = 3
gr.AICc = apply(gr.loglik, 2, function(x) -2*sum(x) + 2*n.df*n.loc + 2*(n.df*n.loc)*(n.df*n.loc+1)/(n.time*n.loc)-n.df*n.loc-1)

## Parameter figures
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          gwgr.slope,
          #interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (mm/month/10 years) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          col=cmap(200),
          zlim = c(-0.2,0.2),
          size = c(2,2), lratio = 0.25,
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ks_elevation$longitude,ks_elevation$latitude,
          gr.slope,
          #interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (mm/month/10 years) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          col=cmap(200),
          zlim = c(-0.1,0.1),
          size = c(2,2), lratio = 0.25,
          legend = "vertical")

gr_z_score = gr.slope/gr.slope.se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ks_elevation$longitude,ks_elevation$latitude,
          gr_z_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z-scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=2),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=cmap(100),
          zlim = c(-2.5,2.5),
          legend = "vertical")

gwgr_z_score = gwgr.slope/gwgr.slope.se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          gwgr_z_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z-scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=2),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=cmap(100),
          zlim = c(-7,7),
          legend = "vertical")

## Supplementary material
### Model selection via 10-fold CV
d_bis = seq(1.5,5,by = 0.5)
set.seed(2021)
n.fold = 10
shuffle.ind = sample(dim(pr_exp)[3])
cv.folds = cut(seq(1,dim(pr_exp)[3]),breaks=n.fold,labels=FALSE)
pred.resi = array(0, dim=c(n.time, n.season, n.loc))
pred.err = matrix(0, nrow =  length(d_bis), ncol = n.fold)
set.seed(2021)
for(nf in 1:n.fold)
{
  print(nf)
  test.ind = shuffle.ind[which(cv.folds==nf)]
  train.ind = shuffle.ind[-which(cv.folds==nf)]
  train.ind = train.ind[order(train.ind, decreasing = FALSE)]
  for(j in 1:length(d_bis))
  {
    w1 = gwr.lm(distance,d_bis[j])
    w1[test.ind,test.ind] = 0
    for(i in 1:n.loc)
    {
      w_i = w1[i,]
      for(t in 1:n.season)
      {
        y = y.all[[t]]
        m_mle2 = bbmle::mle2(gwr_gamma,optimfun = "BFGS",
                             start = list(b0 = 2,b1 = slope_est[i,t], b2=0.01,
                                          logshape = log(shape_est[i,t])))
        beta.t = coef(m_mle2)[1:3]
        pred.resi[,t,i] = exp(x_mult[,,i]%*%beta.t)[,1]
      }
    }
    pred.err[j,nf] = mean(apply(pred.resi[,,test.ind] - pr_exp[,,test.ind], 3, function(x) sum(x^2)))
  }
}


## Model selection via AIC
d_bis = seq(1,10,by = 0.5)
gwgr.AICc.lm = matrix(0, nrow = length(d_bis), ncol = n.season)
gwgr.loglik.t = matrix(0, nrow = n.loc, ncol = n.season)
for(j in 1:length(d_bis))
{
  print(j)
  w1 = gwr.lm(distance,d_bis[j])
  for(i in 1:n.loc)
  {
    w_i = w1[i,]
    for(t in 1:n.season)
    {
      y = y.all[[t]]
      m_mle2 = bbmle::mle2(gwr_gamma,optimfun = "BFGS",
                           start = list(b0 = apply(gr.b0,2,mean)[t],b1 = slope_est[i,t], b2=apply(gr.h,2,mean)[t],
                                        logshape = log(shape_est[i,t])))
      gwgr.loglik.t[i,t] = gwr_loglik(coef(m_mle2)[1],coef(m_mle2)[2],coef(m_mle2)[3],coef(m_mle2)[4])
    }
  }
  gwgr.AICc.lm[j,] = apply(gwgr.loglik.t, 2, function(x) -2*sum(x) + 2*n.df*n.loc + 2*(n.df*n.loc)*(n.df*n.loc+1)/(n.time*n.loc)-n.df*n.loc-1)
}#

