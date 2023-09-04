###########################################################################################
## This file provides the code for producing suppelmentary figures/tables included in the manuscript
###########################################################################################

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
library(car)
## color map for parameters
#cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
cmap = colorRampPalette(rev((brewer.pal(11, "RdBu"))[3:9])) 
rev.cmap = colorRampPalette(brewer.pal(11, "RdBu")[3:9])
neg_cmap = colorRampPalette(rev(brewer.pal(11, "RdBu"))[3:6]) ## all blue
pos_cmap = colorRampPalette(rev((brewer.pal(11, "RdBu"))[3:6])) ## all red
## color map for data
pr.cmap <- colorRampPalette(brewer.pal(11, "Spectral")[6:11])
tas.cmap <- colorRampPalette(rev(brewer.pal(11, "Spectral"))[1:6])
ele.cmap <- colorRampPalette(terrain.colors(10, alpha = 1))

code_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/R/"
data_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Data" 
res_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Results"
pic_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Figures"
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
load(paste0(res_path,"ca_data.rdata"))
load(paste0(res_path,"ca_slope_v_new.rdata"))
load(paste0(res_path,"ca_pre_coef.rdata"))
load(paste0(res_path,"ca_tas_coef.rdata"))
ca_elevation = ca_data$ca_elevation
season.var = c("winter","spring","summer","fall" )
n.season = length(season.var)

plot_list = function(x, main = "", xlab = ""){
  maxy = sapply(x, function(j) j$y)
  plot(x[[1]], ylim = range(maxy), xlab = xlab, main = main,# xlim = c(-3,3),
       col = "darkgray", lwd = 0.5)
  for (j in seq_along(x)) {
    lines(x[[j]], col = "darkgray", lwd = 0.5)
  }
}

##########################################
#### Supplementary Figures for CA ####
##########################################
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

time_dat = data.frame(year = rep(1950:2005,each=12),
                      month = rep(1:12,56),
                      ts.point = 1:(12*56))  ## historical measures
time_rcp = data.frame(year = rep(2006:2100,each=12),
                      month = rep(1:12,95),
                      ts.point = 1:(12*95))
## temperature
gls.tas.summary = ca_data$gls.tas.summary
gls.tas.rcp = ca_data$gls.tas.rcp

## precipitation
gls.pr.summary = ca_data$gls.pr.summary
gls.pr.rcp = ca_data$gls.pr.rcp

n.loc = length(unique(gls.pr.summary$loc))
time.yr = (range(gls.pr.summary$year)[1] : range(gls.pr.summary$year)[2])
time.rcp.yr = (range(gls.pr.rcp$year)[1] : range(gls.pr.rcp$year)[2])
#time.pts = c(time.yr - mean(time.yr))/sd(time.yr)
time.pts = c(time.yr - mean(time.yr))/10
time.rcp.pts = c(time.rcp.yr - mean(time.rcp.yr))/10
n.time = length(time.pts)
n.rcp.time = length(time.rcp.pts)

gls.tas.summary$time.pts = time.pts
gls.tas.rcp$time.pts = time.rcp.pts
gls.tas.all = rbind(gls.tas.summary,gls.tas.summary)

gls.pr.summary$time.pts = time.pts
gls.pr.rcp$time.pts = time.rcp.pts
gls.pr.all = rbind(gls.pr.summary,gls.pr.rcp)


### 3. Correlation checking
##### Temperature
n.lon = length(unique(gls.tas.summary$lon))
n.lat = length(unique(gls.tas.summary$lat))
simgrid = cbind(as.numeric(as.factor(ca_elevation$longitude)),
                as.numeric(as.factor(ca_elevation$latitude)))
distance = as.matrix(dist(simgrid))
r.distance = as.matrix(dist(simgrid)) * 50 ##0.44 degree = 50km

## correlation linear/quadratic test
resids = ca_tas_coef$ca_tas_nw$resids
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
  
  plot(tmp.dat$vec.dist, tmp.dat$tmp.cor,pch=".", cex=0.8,
       main = season.var[t], xlab = "distance",ylab = "correlation")
  pred.cor = predict(fit_lm1,newdata = data.frame(vec.dist = dist.seq))
  lines(dist.seq, pred.cor+1, col="red")
  pred.cor = predict(fit_lm21,newdata = data.frame(dist2 = dist.seq^2))
  #lines(dist.seq, (1-pred.cor^2), col="blue")
  lines(dist.seq, (1-pred.cor)^2, col="blue")
  legend("topright",c("Linear","Bisqaure"), col=c("red","blue"),lty=1, cex = 0.6)
  
  plot(tmp.dat$vec.dist, tmp.dat$log.cor,pch=".", cex=0.8,
       main = season.var[t], xlab = "distance",ylab = "log(cor)")
  abline(fit_loglm, col="red")
  pred.cor = predict(fit_logquad,newdata = data.frame(dist2 = dist.seq^2))
  lines(dist.seq, pred.cor, col="blue")
  legend("topright",c("Exponential","Gaussian"), col=c("red","blue"),lty=1, cex = 0.6)
  
  #tmp.quad = lm(tmp.cor ~ vec.dist + I(vec.dist^2), data = tmp.dat)
  #dist = seq(0, max(vec.dist), by=0.1)
  #pred.cor = predict(tmp.quad,newdata = data.frame(vec.dist = dist))
  #lines(dist, pred.cor, col="blue")
  #legend("topright", legend = c("linear","quadratic"),col = c("red","blue"), lty=1)
}

rownames(cor.rsquare) = m.type
colnames(cor.rsquare) = season.var
cor.rsquare

## cross-validation for selecting threshold
load(paste0(res_path,"ca_tas_cv.rdata"))
d_bis = seq(1,10,by = 0.5)

par(mar=c(5, 7, 4, 6) + 0.1)
plot(d_bis[2:11], apply(pred.err,1,mean)[2:11], pch=16, axes=FALSE, 
     ylim=range(apply(pred.err,1,mean)[2:11]), xlab="", ylab="", 
     type="b",col="black", main="")
axis(2, ylim=range(apply(pred.err,1,mean)[2:11]),col="black",las=1)  ## las=1 makes horizontal labels
mtext("Cross-validation MSE",side=2,line=4)
axis(1,pretty(range(d_bis[2:10]),10))
mtext("Correlation extent",side=1,col="black",line=2.5)  

## Add Legend
legend("topleft",legend=c("AIC","Cross-validation"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

## AIC values
## 63776.39, 66233.08 (no elevation), 65921.86 (no cross-sea)

# CoPE temperature
## Normality check
norm_gls_test = norm_gwr_test = matrix(NA, nrow = n.loc, ncol = n.season)
for(s in 1:n.season)
{
  for(i in 1:n.loc)
  {
    norm_gls_test[i,s] = shapiro.test(ca_tas_coef$ca_tas_nw$resids[i,,s])$p.value
    norm_gwr_test[i,s] = shapiro.test(ca_tas_coef$ca_tas_w$residuals[,s,i])$p.value
  }
  #plot_list(d.resids, main = season.var[j])
}

par(mfrow = c(2,2), oma=c(1, 1, 1,1),  mar = c(2, 2, 2, 2))
ecdf.prop = seq(0.01, 0.99, len = 99)
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(norm_gwr_test[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
}

### Autocorrelation checking
lbtests_gls_pval = lbtests_gwr_pval = matrix(NA, nrow = n.loc, ncol = n.season)
for(s in 1:n.season)
{
  for(i in 1:n.loc)
  {
    lbtests_gls_pval[i,s] = Box.test(ca_tas_coef$ca_tas_nw$resids[i,,s],type = "Ljung-Box",lag = 10)$p.value
    lbtests_gwr_pval[i,s] = Box.test(ca_tas_coef$ca_tas_w$residuals[,s,i],type = "Ljung-Box",lag = 10)$p.value
  }
}

par(mfrow = c(2,2), oma=c(1, 1, 1,1),  mar = c(2, 2, 2, 2))
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(lbtests_gwr_pval[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  if(s==4)
    legend("bottomright", bty='n', xpd=NA,cex = 1,
           legend=c("ECDF","Uniform","95% CI"),lty = 1, 
           col=c("black","blue","red"))
}
# legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
#        legend=c("ECDF","Uniform","95% CI"),lty = 1, 
#        col=c("black","blue","red"))

# CoPE precipitation
load("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/ca_pre_dw_p.rdata")
load("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/ca_pre_dresid.rdata")

## Normality check
ecdf.prop = seq(0.01, 0.99, len = 99)
par(mfrow = c(2,2), oma=c(1, 1, 1,1),  mar = c(2, 2, 2, 2))
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(dw_p[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],
       xlab = "x", ylab = "ECDF",type="l",ylim = c(0,1))
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  if(s==4)
    legend("bottomright", bty='n', xpd=NA,cex = 1,
           legend=c("ECDF","Uniform","95% CI"),lty = 1, 
           col=c("black","blue","red"))
}

# legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
#        legend=c("ECDF","Uniform","95% CI"),lty = 1, 
#        col=c("black","blue","red"))

## Autocorrelation check
lbtests_gwr_pval = matrix(NA, nrow = n.loc, ncol = n.season)
for(s in 1:n.season)
{
  for(i in 1:n.loc)
  {
    lbtests_gwr_pval[i,s] = Box.test(gwgr.pre.ca$res.gwgr$dat_resid[,s,i],type = "Ljung-Box",lag = 10)$p.value
  }
}

par(mfrow = c(2,2), oma=c(1, 1, 1,1),  mar = c(2, 2, 2, 2))
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(lbtests_gwr_pval[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],
       xlab = "x", ylab = "ECDF",type="l",ylim = c(0,1))
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  if(s==4)
    legend("bottomright", bty='n', xpd=NA,cex = 1,
           legend=c("ECDF","Uniform","95% CI"),lty = 1, 
           col=c("black","blue","red"))
}
# legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
#        legend=c("ECDF","Uniform","95% CI"),lty = 1, 
#        col=c("black","blue","red"))
