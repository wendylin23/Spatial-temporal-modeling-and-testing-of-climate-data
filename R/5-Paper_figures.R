###########################################################################################
## This file provides the code for producing all figures included in the manuscript
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

code_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/R/"
data_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Data" 
res_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Data/Results/"
pic_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Figures"
setwd(code_path)
source(paste0(code_path,"gwr.R"))
source(paste0(code_path,"varx_fixed.R"))
source(paste0(code_path,"util.R"))

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

## get aspect ratio from lon/lat
map_aspect = function(x, y) {
  x.center = sum(range(x)) / 2
  y.center = sum(range(y)) / 2
  x.dist = ggplot2:::dist_central_angle(x.center + c(-0.5, 0.5), rep(y.center, 2))
  y.dist = ggplot2:::dist_central_angle(rep(x.center, 2), y.center + c(-0.5, 0.5))
  ratio = y.dist / x.dist
  diff(range(y)) / diff(range(x)) * ratio
}

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
load(paste0(res_path,"ca_pre_dw_p.rdata"))
ca_elevation = ca_data$ca_elevation
season.var = c("winter","spring","summer","fall" )
n.season = length(season.var)
loc.names = c("San Francisco","San Diego","Yosemite","Death Valley")
select.loc = data.frame(sf = c(37.75,-122.25), sd = c(32.75,-117.25),
                        yo = c(37.75,-119.25), dv = c(36.25,-116.75))
select.ind = c(126, 1, 132, 87)
select.loc.col = c("red","green","blue","darkgrey")
season.var = c("winter","spring","summer","fall")
sub.loc = list(x = c(-122.25, -117.25, -119.25, -116.75), 
               y = c(37.75,32.75,37.75,36.25),
               labels = loc.names)

### Examples of final Cope results
## Hist - temp
#jpeg(paste0(pic_path,"/Plot3.jpeg"), width = 600, height = 600, units = 'mm', res = 300)
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
Simul_Cope(map = ca_elevation, par = slope_new$slope_hist$beta_w_tau,
           par.se =  slope_new$slope_hist$sigma_w_tau,
           resids = ca_tas_coef$ca_tas_w$residuals,
           level.set = c(0.1,0.15,0.2), zlim = c(0,0.6),
           col.map = pos_cmap(200),
           sub.loc=sub.loc, select.ind = select.ind)
#dev.off()

## Hist - pre
Simul_Cope(map = ca_elevation, par = gwgr.pre.ca$res.gwgr$gwgr.slope,
           par.se =  gwgr.pre.ca$res.gwgr$gwgr.slope.se,
           resids = gwgr.pre.ca$res.gwgr$dat_resid,
           level.set = c(-0.05,0,0.05), zlim = c(-0.2,0.2),
           col.map = rev.cmap(200),type = "prep",
           sub.loc=sub.loc, select.ind = select.ind)

## RCP - temp
Simul_Cope(map = ca_elevation, par = slope_new$slope_rcp$beta_w_rcp_tau,
           par.se =  slope_new$slope_rcp$sigma_w_rcp_tau,
           resids = ca_tas_coef$ca_tas_pre$residuals,
           level.set = c(0.5,0.55,0.6), zlim = c(0.4,1),
           col.map = pos_cmap(200),
           sub.loc=sub.loc, select.ind = select.ind)

## RCP - pre
Simul_Cope(map = ca_elevation, par = gwgr.pre.ca$future.gwgr$gwgr.slope,
           par.se =  gwgr.pre.ca$future.gwgr$gwgr.slope.se,
           resids = gwgr.pre.ca$future.gwgr$dat_resid,
           level.set = c(-0.05,0,0.05), zlim = c(-0.25,0.25),
           col.map = rev.cmap(200),type = "prep",
           sub.loc=sub.loc, select.ind = select.ind)


# map of elevation 1.13:1 (1.2:1, 620*600)
par(oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          ca_elevation$elevation_geonames,
          interp.args = list(no.X = 200, no.Y = 200),
          common.legend = FALSE,
          #proj="albers",par=c(30,40),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = c("Elevation (m) in CA"),
          points = sub.loc,
          points.args = list(pch = 20, col = "blue"),
          text = sub.loc,
          text.args = list(pos = 3, col = "blue",cex=0.9),
          mtext.args = list(cex = 1),
          col=ele.cmap(200),
          legend.axis.args = list(cex.axis=1),
          lratio = 0.2,
          legend = "vertical")

### 2. Data visualizationhttps://www.google.com/search?q=california+map&rlz=1C5CHFA_enUS932US932&sxsrf=APq-WBuSjejIxM20Xqu0eCUyqBy-6a0snQ%3A1649300342678&ei=dlNOYoqKKcyblwSl4K7wAg&ved=0ahUKEwiKnbbI-oD3AhXMzYUKHSWwCy4Q4dUDCA4&uact=5&oq=california+map&gs_lcp=Cgdnd3Mtd2l6EAMyBAgjECcyBAgjECcyBAgjECcyBggAEAcQHjIFCAAQgAQyBggAEAcQHjIGCAAQBxAeMgYIABAHEB4yBggAEAcQHjIGCAAQBxAeOgQIABBDOggIABAHEAoQHjoFCC4QgAQ6CwgAEIAEELEDEIMBOhAILhCxAxCDARDHARDRAxBDOg0IABCABBCHAhCxAxAUOgcIIxCxAhAnOgIIJkoECEEYAEoECEYYAFAAWIYRYLASaABwAXgAgAGiAYgBjgySAQQwLjEymAEAoAEBwAEB&sclient=gws-wiz#
#### Temperature
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
          text.args = list(pos = 3, col = "darkred",cex=1.1),
          zlim = c(-12,32),
          col=tas.cmap(200),
          size = c(2, 2), lratio = 0.2,
          legend = "vertical")

## temporal seasonal trend (plus future)
par(mfrow = c(2,2), oma=c(1, 1, 2,1),  mar = c(2, 2, 2, 2))
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
  if(s==2)
    legend("bottomright", bty='n', xpd=NA,cex = 0.9,
           legend=season.var,lty = 1, col=1:length(season.var))
}
#legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
#       legend=season.var,lty = 1, col=1:length(season.var))
mtext("Average seasonal temperature (C) at four typical locations", outer = TRUE, cex = 1.2)

#### Precipitation
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


### 3. Temporal model checking
##### Temperature

#### 3.1.1 Normality check
plot_list = function(x, main = "", xlab = ""){
  maxy = sapply(x, function(j) j$y)
  plot(x[[1]], ylim = range(maxy), xlab = xlab, main = main,# xlim = c(-3,3),
       col = "darkgray", lwd = 0.5)
  for (j in seq_along(x)) {
    lines(x[[j]], col = "darkgray", lwd = 0.5)
  }
}

norm_test = matrix(0, nrow = n.loc, ncol = n.season)

for(j in 1:n.season)
{
  d.resids = list()
  for(i in 1:n.loc)
  {
    d.resids[[i]] = density(gls.tas.summary[gls.tas.summary$loc==i,season.var[j]])
    norm_test[i,j] = shapiro.test(gls.tas.summary[gls.tas.summary$loc==i,season.var[j]])$p.value
  }
  #plot_list(d.resids, main = season.var[j])
}

normal.res = table(norm_test < 0.05)
names(normal.res) = c("normal","non-normal")
normal.res

par(mfrow = c(2,2), oma=c(1, 1, 1,5),  mar = c(2, 2, 2, 2))
ecdf.prop = seq(0.01, 0.99, len = 99)
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(norm_test[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c("ECDF","Uniform","95% CI"),lty = 1, 
       col=c("black","blue","red"))
mtext("Empirical CDFs of the p values from Shapiro-Wilk test of normality", outer = TRUE, cex = 1)

#### 3.1.2 Autocorrelation check
lm_acf = lm_pacf = matrix(NA, nrow = n.loc, ncol = n.season)
resid_mat = array(0, dim = c(n.loc, n.season, n.time))
acf_sig = qnorm(1 - 0.05 / 2) / sqrt(n.time)
box_pval = matrix(0, nrow = n.loc, ncol = n.season)

for(i in 1:n.loc)
{
  for(j in 1:n.season)
  {
    tas_exp = gls.tas.summary[gls.tas.summary$loc==i,season.var[j]]
    #tas_out = tas_exp[-1]
    #tas_lag1 = tas_exp[-n.time]
    lm.fit = lm(tas_exp ~ time.pts)
    resid_mat[i,j,] = residuals(lm.fit)
    acf.res = acf(residuals(lm.fit),plot = FALSE)
    lm_acf[i,j] = acf.res$acf[1,,]
    pacf.res = pacf(residuals(lm.fit),plot = FALSE)
    lm_pacf[i,j] = pacf.res$acf[1,,]
    box_pval[i,j] = Box.test(resid_mat[i,j,],type = "Ljung-Box",lag = 10)$p.value
  }
}

par(mfrow = c(1,1))
boxplot(lm_pacf, ylim = c(-acf_sig-0.1,acf_sig+0.1),
        names=season.var,ylab = "autocorrelation")
abline(h = acf_sig,col="red",lty=2)
abline(h = -acf_sig,col="red",lty=2)

par(mfrow = c(2,2), oma=c(1, 1, 1,5),  mar = c(2, 2, 2, 2))
ecdf.prop = seq(0.01, 0.99, len = 99)
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(box_pval[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c("ECDF","Uniform","95% CI"),lty = 1, 
       col=c("black","blue","red"))
mtext("Empirical CDFs of the p values from Ljung-Box test of autocorrelation", outer = TRUE, cex = 1)

#### 3.1.3 Cross-season checking
multi.coef = multi.pval = matrix(0, nrow = n.loc, ncol = n.season)
for(i in 1:n.loc)
{
  tas_exp = gls.tas.summary[gls.tas.summary$loc==i,season.var]
  # Step 1: seasonal correlation
  season.sub.m1 = lm(tas_exp$spring ~ tas_exp$winter)
  season.sub.m2 = lm(tas_exp$summer ~ tas_exp$spring)
  season.sub.m3 = lm(tas_exp$fall ~ tas_exp$summer)
  season.sub.m4 = lm(tas_exp$winter[2:n.time] ~ tas_exp$fall[1:(n.time-1)])
  multi.coef[i,] = round(c(coef(season.sub.m1)[2],coef(season.sub.m2)[2],
                           coef(season.sub.m3)[2],coef(season.sub.m4)[2]),4)
  multi.pval[i,] = round(c(summary(season.sub.m1)$coefficients[2,4],
                           summary(season.sub.m2)$coefficients[2,4],
                           summary(season.sub.m3)$coefficients[2,4],
                           summary(season.sub.m4)$coefficients[2,4]),4)
}
colnames(multi.coef) = colnames(multi.pval) = c("sp ~ wi","su ~ sp","fa ~ su","wi ~ fa")

par(mfrow=c(1,1))
boxplot(multi.pval,ylab = "p-value of coefficients")
abline(h=0.05,col="red")

par(mfrow = c(2,2), oma=c(1, 1, 1,5),  mar = c(2, 2, 2, 2))
ecdf.prop = seq(0.01, 0.99, len = 99)
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(multi.pval[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],type="l", ylim=c(0,1))
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c("ECDF","Uniform","95% CI"),lty = 1, 
       col=c("black","blue","red"))
mtext("Empirical CDFs of the p values from Ljung-Box test of autocorrelation", outer = TRUE, cex = 1)

#### Precipitation
#### 3.2.2 Autocorrelation check
par(mfrow = c(2,2), oma=c(1, 1, 1,5),  mar = c(2, 2, 2, 2))
ecdf.prop = seq(0.01, 0.99, len = 99)
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(dw_p[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c("ECDF","Uniform","95% CI"),lty = 1, 
       col=c("black","blue","red"))
mtext("Empirical CDFs of the p values from Durbin Watson test of autocorrelation", outer = TRUE, cex = 1)

#### 3.2.3 cross-seasonal effect
multi.coef = multi.pval = matrix(0, nrow = n.loc, ncol = n.season)
for(i in 1:n.loc)
{
  tas_exp = gls.pr.summary[gls.pr.summary$loc==i,season.var]
  # Step 1: seasonal correlation
  season.sub.m1 = lm(tas_exp$spring ~ tas_exp$winter)
  season.sub.m2 = lm(tas_exp$summer ~ tas_exp$spring)
  season.sub.m3 = lm(tas_exp$fall ~ tas_exp$summer)
  season.sub.m4 = lm(tas_exp$winter[2:n.time] ~ tas_exp$fall[1:(n.time-1)])
  multi.coef[i,] = round(c(coef(season.sub.m1)[2],coef(season.sub.m2)[2],
                           coef(season.sub.m3)[2],coef(season.sub.m4)[2]),4)
  multi.pval[i,] = round(c(summary(season.sub.m1)$coefficients[2,4],
                           summary(season.sub.m2)$coefficients[2,4],
                           summary(season.sub.m3)$coefficients[2,4],
                           summary(season.sub.m4)$coefficients[2,4]),4)
}
colnames(multi.coef) = colnames(multi.pval) = c("sp ~ wi","su ~ sp","fa ~ su","wi ~ fa")

par(mfrow=c(1,1))
boxplot(multi.pval,ylab = "p-value of coefficients")
abline(h=0.05,col="red")

#### 4 temporal modeling

### temperature (time trend variable)
beta_est = ca_tas_coef$ca_tas_nw$beta_est
beta_se = ca_tas_coef$ca_tas_nw$beta_se
z_tas_score = beta_est/beta_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          beta_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes in deterministic time trend (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=cmap(200),
          zlim = c(-0.5,0.5),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          z_tas_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of slopes in deterministic time trend of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.2,
          col=cmap(200),
          zlim = c(-12,12),
          legend = "vertical")

## transformed slopes
slope_est = slope_new$slope_hist$beta_tau
slope_se = slope_new$slope_hist$sigma_tau
z_slope_score = slope_est/slope_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          slope_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes in regression time trend (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,0.6),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          z_slope_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of slopes in regression time trend of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,6),
          legend = "vertical")


### precipitation
slope_est = gwgr.pre.ca$res.gr$gwgr.slope
slope_se = gwgr.pre.ca$res.gr$gwgr.slope.se
z_pre_score = slope_est/slope_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          slope_est,
          #interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (mm/month/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          col=rev.cmap(200),
          zlim = c(-0.2,0.2),
          size = c(2,2), lratio = 0.2,
          legend = "vertical")

z_pre_score = ifelse(!is.na(z_pre_score),z_pre_score,0)
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          z_pre_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z-scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=rev.cmap(200),
          zlim = c(-7,7),
          legend = "vertical")

#### 5 spatial-temporal modeling
### temperature
w_beta_est = t(ca_tas_coef$ca_tas_w$coef[2,,])
w_beta_se = t(ca_tas_coef$ca_tas_w$se.coef[2,,])
zw_tas_score = w_beta_est/w_beta_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          w_beta_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes in deterministic time trend (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=cmap(200),
          zlim = c(-0.5,0.5),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          zw_tas_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of slopes in deterministic time trend of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.2,
          col=cmap(200),
          zlim = c(-12,12),
          legend = "vertical")

w_slope_est = slope_new$slope_hist$beta_w_tau
w_slope_se = slope_new$slope_hist$sigma_w_tau
zw_slope_score = w_slope_est/w_slope_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          w_slope_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes in regression time trend (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,0.6),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          zw_slope_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of slopes in deterministic time trend of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,12),
          legend = "vertical")

## weighted vs. non-weighted model
par(mfrow = c(2,2), oma=c(1, 1, 2,1),  mar = c(2, 2, 2, 2))
col1 = rgb(1,0,0,1/4)
col2 = rgb(0,1,0,1/4)
for(i in 1:n.season)
{
  tmp.z = cbind(z_tas_score[,i],zw_tas_score[,i])
  hist(tmp.z[,1],xlab = "Z-score",main = season.var[i],
       seq(from=-20, to=20, by=0.5),probability = TRUE,
       col=col1,xlim = c(-15,15))
  #hist(z_score$slope_m1[,2],col=rgb(0,1,1,1/4), add = T, breaks =30)
  hist(tmp.z[,2],col=col2,probability = TRUE, add = T,seq(from=-20, to=20, by=0.5))
  if(i==n.season)
    legend("topright",legend = c("MTSR","GWMTSR"),
           fill = c(col1,col2),
           bty = 'n',cex = 0.75,
           border = NA)
}
mtext("Histogram of slope z-scores of deterministic time trend", outer = TRUE, cex = 1.2)

par(mfrow = c(2,2), oma=c(1, 1, 2,1),  mar = c(2, 2, 2, 2))
col1 = rgb(1,0,0,1/4)
col2 = rgb(0,1,0,1/4)
for(i in 1:n.season)
{
  tmp.z = cbind(z_slope_score[,i],zw_slope_score[,i])
  hist(tmp.z[,1],xlab = "Z-score",main = season.var[i],
       seq(from=-20, to=20, by=0.5),probability = TRUE,
       col=col1,xlim = c(-15,15))
  #hist(z_score$slope_m1[,2],col=rgb(0,1,1,1/4), add = T, breaks =30)
  hist(tmp.z[,2],col=col2,probability = TRUE, add = T,seq(from=-20, to=20, by=0.5))
  if(i==n.season)
    legend("topright",legend = c("MTSR","GWMTSR"),
           fill = c(col1,col2),
           bty = 'n',cex = 0.75,
           border = NA)
}
mtext("Histogram of slope z-scores of regression time trend", outer = TRUE, cex = 1.2)


### precipitation
gwgr.slope = gwgr.pre.ca$res.gwgr$gwgr.slope
gwgr.slope.se = gwgr.pre.ca$res.gwgr$gwgr.slope.se
zw_pre_score = gwgr.slope/gwgr.slope.se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          gwgr.slope,
          #interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (mm/month/10 years) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          col=rev.cmap(200),
          zlim = c(-0.2,0.2),
          size = c(2,2), lratio = 0.2,
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          zw_pre_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z-scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=rev.cmap(200),
          zlim = c(-7,7),
          legend = "vertical")

par(mfrow = c(2,2), oma=c(1, 1, 2,1),  mar = c(2, 2, 2, 2))
col1 = rgb(1,0,0,1/4)
col2 = rgb(0,1,0,1/4)
for(i in 1:n.season)
{
  tmp.z = cbind(z_pre_score[,i],zw_pre_score[,i])
  hist(tmp.z[,1],xlab = "Z-score",main = season.var[i],
       seq(from=-20, to=20, by=0.5),probability = TRUE, 
       col=col1,xlim = c(-8,8))
  #hist(z_score$slope_m1[,2],col=rgb(0,1,1,1/4), add = T, breaks =30)
  hist(tmp.z[,2],col=col2, add = T,
       seq(from=-20, to=20, by=0.5), probability = TRUE)
  if(i==n.season)
    legend("topright",legend = c("GR","GWGR"),
           fill = c(col1,col2),
           bty = 'n',cex = 0.75,
           border = NA)
}
mtext("Histogram of slope z-scores", outer = TRUE, cex = 1.2)

### future temperature
fut_beta_est = t(ca_tas_coef$ca_tas_pre$coef[2,,])
fut_beta_se = t(ca_tas_coef$ca_tas_pre$se.coef[2,,])
zw_fut_score = fut_beta_est/fut_beta_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          fut_beta_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,1.1),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          zw_fut_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.2,
          col=cmap(200),
          zlim = c(-12,12),
          legend = "vertical")

fut_slope_est = slope_new$slope_rcp$beta_w_rcp_tau
fut_slope_se = slope_new$slope_rcp$sigma_w_rcp_tau
zw_fut_slope_score = fut_slope_est/fut_slope_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          fut_slope_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0.5,1),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          zw_fut_slope_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          points = sub.loc,
          points.args = list(pch = 16, col = "black",cex=1.1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,50),
          legend = "vertical")

#### 6 CoPE inference
### temperature
### 6.1 Continuity checking
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          ca_tas_coef$ca_tas_nw$beta_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=cmap(200),
          zlim = c(-0.5,0.5),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          t(ca_tas_coef$ca_tas_w$coef[2,,]),
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=cmap(200),
          zlim = c(-0.5,0.5),
          legend = "vertical")

## transformed slopes
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          slope_new$slope_hist$beta_tau,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,0.6),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          slope_new$slope_hist$beta_w_tau,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,0.6),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          slope_new$slope_hist$beta_lm,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.2,
          col=pos_cmap(200),
          zlim = c(0,0.6),
          legend = "vertical")

#### Cope set for CA
#Confidence level for CoPE sets.
alpha_CoPE = 0.1
#Nominal expected false area ratio for FARE sets.
alpha_FARE = 0.05
#Number of realizations to produce in the MC simulations.
N=1000
lon = unique(ca_elevation$longitude)
lon = lon[order(lon)]
lat = unique(ca_elevation$latitude)
x.grid = seq(0,1, length = length(lon))
y.grid = seq(0,1, length = length(lat))
n = dim(ca_tas_coef$ca_tas_w$residuals)[1]

### 6.2 Normality checking
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

par(mfrow = c(2,2), oma=c(1, 1, 1,5),  mar = c(2, 2, 2, 2))
ecdf.prop = seq(0.01, 0.99, len = 99)
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(norm_gls_test[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c("ECDF","Uniform","95% CI"),lty = 1, 
       col=c("black","blue","red"))

par(mfrow = c(2,2), oma=c(1, 1, 1,5),  mar = c(2, 2, 2, 2))
ecdf.prop = seq(0.01, 0.99, len = 99)
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(norm_gwr_test[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c("ECDF","Uniform","95% CI"),lty = 1, 
       col=c("black","blue","red"))

### 6.3 Autocorrelation checking
lbtests_gls_pval = lbtests_gwr_pval = matrix(NA, nrow = n.loc, ncol = n.season)
for(s in 1:n.season)
{
  for(i in 1:n.loc)
  {
    lbtests_gls_pval[i,s] = Box.test(ca_tas_coef$ca_tas_nw$resids[i,,s],type = "Ljung-Box",lag = 10)$p.value
    lbtests_gwr_pval[i,s] = Box.test(ca_tas_coef$ca_tas_w$residuals[,s,i],type = "Ljung-Box",lag = 10)$p.value
  }
}

par(mfrow = c(2,2), oma=c(1, 1, 1,5),  mar = c(2, 2, 2, 2))
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(lbtests_gls_pval[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c("ECDF","Uniform","95% CI"),lty = 1, 
       col=c("black","blue","red"))

par(mfrow = c(2,2), oma=c(1, 1, 1,5),  mar = c(2, 2, 2, 2))
for(s in 1:n.season)
{
  ecdfs = sapply(ecdf.prop, function(y) mean(lbtests_gwr_pval[,s] <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = season.var[s],type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
}
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=c("ECDF","Uniform","95% CI"),lty = 1, 
       col=c("black","blue","red"))
