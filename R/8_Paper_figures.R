###### Packages and functions #####
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
source("~/Dropbox/UCSD/Thesis/3.Precipitation/Code/util.R")
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

res_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/summary_results_20211101/"
code_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/"
pic_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Figures/20220315_supp"
season.var = c("winter","spring","summer","fall" )
n.season = length(season.var)


## get aspect ratio from lon/lat
map_aspect = function(x, y) {
  x.center = sum(range(x)) / 2
  y.center = sum(range(y)) / 2
  x.dist = ggplot2:::dist_central_angle(x.center + c(-0.5, 0.5), rep(y.center, 2))
  y.dist = ggplot2:::dist_central_angle(rep(x.center, 2), y.center + c(-0.5, 0.5))
  ratio = y.dist / x.dist
  diff(range(y)) / diff(range(x)) * ratio
}

### 1. Elevation map of three states
load(paste0(res_path,"ks_data.rdata"))
ks_elevation = ks_data$ks_elevation
season.var = c("winter","spring","summer","fall" )
loc.names = c("Kansas City","Wichita","Oakley")
select.loc = data.frame(kc = c(39.25,-94.75), wi = c(37.75,-97.25),
                        oa = c(39.25,-100.75))
select.ind = c(75,25,63)
select.loc.col = c("red","green","blue")
sub.loc = list(x = c(-94.97,-97.25,-100.75), 
               y = c(39.25,37.75,39.25),
               labels = loc.names)
load(paste0(res_path,"ks_slope_v_new.rdata"))
load(paste0(res_path,"ks_pre_coef.rdata"))
load(paste0(res_path,"ks_tas_coef.rdata"))

## Hist - temp
Simul_Cope(map = ks_elevation, par = slope_new$slope_hist$beta_w_tau,
           par.se =  slope_new$slope_hist$sigma_w_tau,
           resids = ks_tas_coef$ks_tas_w$residuals,
           level.set = c(0.1,0.15,0.2), zlim = c(-0.3,0.3),
           col.map = cmap(200),
           sub.loc=sub.loc, select.ind = select.ind)

## Hist - pre
Simul_Cope(map = ks_elevation, par = gwgr.pre.ks$res.gwgr$gwgr.slope,
           par.se =  gwgr.pre.ks$res.gwgr$gwgr.slope.se,
           resids = gwgr.pre.ks$res.gwgr$dat_resid,
           level.set = c(-0.05,0,0.05), zlim = c(-0.09,0.09),
           col.map = rev.cmap(200),type = "prep",
           sub.loc=sub.loc, select.ind = select.ind)

## RCP - temp
Simul_Cope(map = ks_elevation, par = slope_new$slope_rcp$beta_w_rcp_tau,
           par.se =  slope_new$slope_rcp$sigma_w_rcp_tau,
           resids = ks_tas_coef$ks_tas_pre$residuals,
           level.set = c(0.5,0.55,0.6), zlim = c(0.2,0.8),
           col.map = pos_cmap(200),
           sub.loc=sub.loc, select.ind = select.ind)

## RCP - pre
Simul_Cope(map = ks_elevation, par = gwgr.pre.ks$future.gwgr$gwgr.slope,
           par.se =  gwgr.pre.ks$future.gwgr$gwgr.slope.se,
           resids = gwgr.pre.ks$future.gwgr$dat_resid,
           level.set = c(-0.05,0,0.05), zlim = c(-0.06,0.06),
           col.map = rev.cmap(200),type = "prep",
           sub.loc=sub.loc, select.ind = select.ind)


## 0.46:1 (0.51:1, 600*312)
c_adj = cosd(mean(ks_elevation$latitude))
p_adj = diff(range(ks_elevation$longitude)) * c_adj/diff(range(ks_elevation$latitude))

par(oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))
autoimage(ks_elevation$longitude,ks_elevation$latitude,
          ks_elevation$elevation_geonames,
          interp.args = list(no.X = 200, no.Y = 200),
          common.legend = FALSE,
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = c("Elevation (m) in KS"),
          points = sub.loc,
          points.args = list(pch = 20, col = "blue"),
          text = sub.loc,
          text.args = list(pos = 3, col = "blue",cex=1.1),
          mtext.args = list(cex = 1),
          col=ele.cmap(200),
          lratio = 0.2,
          legend.axis.args = list(cex.axis=1),
          legend = "vertical")

## CO
load(paste0(res_path,"co_data.rdata"))
season.var = c("winter","spring","summer","fall" )
loc.names = c("Denver","Aspen","Grand Junction")
select.loc = data.frame(dv = c(39.75,-104.75), ap = c(39.25,-106.75),
                        gj = c(39.25,-108.75))
select.ind = c(79,61,57)
select.loc.col = c("red","green","blue")
sub.loc = list(x = c(-104.75,-106.75,-108.75), 
               y = c(39.75,39.25,39.25),
               labels = loc.names)
co_elevation = co_data$co_elevation

load(paste0(res_path,"co_slope_v_new.rdata"))
load(paste0(res_path,"co_pre_coef.rdata"))
load(paste0(res_path,"co_tas_coef.rdata"))

## Hist - temp (600*435)
Simul_Cope(map = co_elevation, par = slope_new$slope_hist$beta_w_tau,
           par.se =  slope_new$slope_hist$sigma_w_tau,
           resids = co_tas_coef$co_tas_w$residuals,
           level.set = c(0.1,0.15,0.2), zlim = c(-0.3,0.3),
           col.map = cmap(200),
           sub.loc=sub.loc, select.ind = select.ind)

## Hist - pre
Simul_Cope(map = co_elevation, par = gwgr.pre.co$res.gwgr$gwgr.slope,
           par.se =  gwgr.pre.co$res.gwgr$gwgr.slope.se,
           resids = gwgr.pre.co$res.gwgr$dat_resid,
           level.set = c(-0.05,0,0.05), zlim = c(-0.05,0.05),
           col.map = rev.cmap(200),type = "prep",
           sub.loc=sub.loc, select.ind = select.ind)

## RCP - temp
Simul_Cope(map = co_elevation, par = slope_new$slope_rcp$beta_w_rcp_tau,
           par.se =  slope_new$slope_rcp$sigma_w_rcp_tau,
           resids = co_tas_coef$co_tas_pre$residuals,
           level.set = c(0.5,0.55,0.6), zlim = c(0.2,0.8),
           col.map = pos_cmap(200),
           sub.loc=sub.loc, select.ind = select.ind)

## RCP - pre
Simul_Cope(map = co_elevation, par = gwgr.pre.co$future.gwgr$gwgr.slope,
           par.se =  gwgr.pre.co$future.gwgr$gwgr.slope.se,
           resids = gwgr.pre.co$future.gwgr$dat_resid,
           level.set = c(-0.05,0,0.05), zlim = c(-0.07,0.07),
           col.map = rev.cmap(200),type = "prep",
           sub.loc=sub.loc, select.ind = select.ind)


##0.69:1 (0.74:1)
par(oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          co_elevation$elevation_geonames,
          interp.args = list(no.X = 200, no.Y = 200),
          common.legend = FALSE,
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = c("Elevation (m) in CO"),
          points = sub.loc,
          points.args = list(pch = 20, col = "blue"),
          text = sub.loc,
          text.args = list(pos = 3, col = "blue",cex=1.1),
          mtext.args = list(cex = 1),
          col=ele.cmap(200),
          legend.axis.args = list(cex.axis=1),
          lratio = 0.2,
          legend = "vertical")

load(paste0(res_path,"ca_data.rdata"))
ca_elevation = ca_data$ca_elevation
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

load(paste0(res_path,"ca_slope_v_new.rdata"))
load(paste0(res_path,"ca_pre_coef.rdata"))
load(paste0(res_path,"ca_tas_coef.rdata"))

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
load(paste0(res_path,"ca_pre_coef.rdata"))
load(paste0(res_path,"ca_tas_coef.rdata"))
#load(paste0(res_path,"ca_tas_cv.rdata"))
#load(paste0(res_path,"ca_tas_aic_w2.rdata"))

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
load("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/ca_pre_dw_p.rdata")

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
load(paste0(res_path,"ca_slope_v_new.rdata"))
load(paste0(res_path,"ca_pre_coef.rdata"))
load(paste0(res_path,"ca_tas_coef.rdata"))

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

### 6.4 CoPE set
#### OLS
pdf(file = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Figures/20220118/Cope.pdf")
level.set = c(0.05,0.1,0.15)
n = dim(ca_tas_coef$ca_tas_w$residuals)[1]
for(l in 1:length(level.set))
{
  par(mfrow = c(2,2), oma = c(0, 0, 0, 3), mar = c(2, 1, 1, 2))
  zlim = c(-0.5,0.5)
  for(s in 1:n.season)
  {
    
    beta.long = data.frame(lon = ca_elevation$longitude, 
                           lat = ca_elevation$latitude,
                           b.gls = ca_tas_coef$ca_tas_nw$beta_est[,s],
                           b.gwr = ca_tas_coef$ca_tas_w$coef[2,s,]
    )
    se.long = data.frame(lon = ca_elevation$longitude, 
                         lat = ca_elevation$latitude,
                         b.gls = ca_tas_coef$ca_tas_nw$beta_se[,s],
                         b.gwr = ca_tas_coef$ca_tas_w$se.coef[2,s,]
    )
    
    beta.gls.wide = acast(beta.long, lon~lat, value.var='b.gls')
    beta.gwr.wide = acast(beta.long, lon~lat, value.var='b.gwr')
    se.gls.wide = acast(se.long, lon~lat, value.var='b.gls')
    se.gwr.wide = acast(se.long, lon~lat, value.var='b.gwr')
    mask = (!is.na(beta.gwr.wide))*1
    
    #Compute the residuals.
    R.gls= R.gwr = array(0,c(length(lon),length(lat),n))
    for(i in 1:n)
    {
      long.t = data.frame(lon = ca_elevation$longitude, 
                          lat = ca_elevation$latitude,
                          #R.gls = ks_tas_coef$ks_tas_nw$resids[,i,2],
                          R.gls = ca_tas_coef$ca_tas_nw$resids[,i,s],
                          R.gwr = ca_tas_coef$ca_tas_w$residuals[i,s,])
      R.gls[,,i] = acast(long.t, lon~lat, value.var='R.gls')
      R.gwr[,,i] = acast(long.t, lon~lat, value.var='R.gwr')
    }
    
    sigma.gls.hat = apply(R.gls,1:2,function(x) sd(x, na.rm = TRUE))
    R.gls.tilde = R.gls / rep(sigma.gls.hat,n)
    sigma.gwr.hat = apply(R.gwr,1:2,function(x) sd(x, na.rm = TRUE))
    R.gwr.tilde = R.gwr / rep(sigma.gwr.hat,n)
    
    sigma.gls.long = apply(ca_tas_coef$ca_tas_nw$resids[,,s],2,function(x) sd(x, na.rm = TRUE))
    R.gls.long = ca_tas_coef$ca_tas_nw$resids[,,s]/sigma.gls.long
    R.gls.sub = t(R.gls.long[select.ind,])
    
    sigma.gwr.long = apply(ca_tas_coef$ca_tas_w$residuals[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gwr.long = ca_tas_coef$ca_tas_w$residuals[,s,]/sigma.gwr.long
    R.gwr.sub = R.gwr.long[,select.ind]
    
    A = (mask==1)
    A[is.na(A)] = FALSE
    for(i in 1:n) R.gls.tilde[,,i] = R.gls.tilde[,,i] * A
    for(i in 1:n) R.gwr.tilde[,,i] = R.gwr.tilde[,,i] * A
    
    #Compute quantile (bootstrap).
    P_MC_gls = MC_gauss(R.gls.tilde,N)
    P_MC_gwr = MC_gauss(R.gwr.tilde,N)
    a_CoPE_gls = a_CoPE_gwr = 0
    while(P_MC_gls(a_CoPE_gls)>alpha_CoPE/3) a_CoPE_gls = a_CoPE_gls+0.01
    while(P_MC_gwr(a_CoPE_gwr)>alpha_CoPE/3) a_CoPE_gwr = a_CoPE_gwr+0.01
    
    #Normalized function.
    level= level.set[l]
    # norm_diff_gls = (beta.gls.wide - level) / se.gls.wide
    # norm_diff_gwr = (beta.gwr.wide - level) / se.gwr.wide
    
    n.lon = length(lon)
    n.lat = length(lat)
    
    norm_diff_gls = (beta.long$b.gls - level) / se.long$b.gls
    norm_diff_gwr = (beta.long$b.gwr - level) / se.long$b.gwr
    
    # image.map(lon,lat,beta.gls.wide,mask=mask,col = tim.colors(64),
    #           horizontal=FALSE,ylab='',xlab='',
    #           main= paste0("GLS ",season.var[s]," level=",level))
    
    slope_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                  ca_tas_coef$ca_tas_nw$beta_est[,s]), 200, 200)$xyz.est
    gls_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                norm_diff_gls), 200, 200)$xyz.est
    gwr_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                norm_diff_gwr), 200, 200)$xyz.est
    mask_large = (!is.na(slope_interp$z))*1
    image.map(slope_interp$x,slope_interp$y,
              slope_interp$z,mask=mask_large,col = cmap(200),
              #horizontal=FALSE,
              legend = FALSE,
              ylab='',xlab='',
              zlim = zlim,
              main= paste0(season.var[s]))
    maps::map("county",add = TRUE)
    points(sub.loc$x, sub.loc$y,pch = 16, col = "black",cex=1.1)
    drawContour(slope_interp$x,slope_interp$y,bound = FALSE,
                slope_interp$z,c=level,col="purple",lty=1)
    drawContour(gls_interp$x,gls_interp$y,gls_interp$z,c=a_CoPE_gls,col="darkred")
    drawContour(gls_interp$x,gls_interp$y,gls_interp$z,c=-a_CoPE_gls,col="darkgreen")
  }
  par(mfrow=c(1, 1), mar=c(1, 1, 0.7, 2), new=FALSE)
  image.plot(zlim= zlim, legend.only=TRUE, col=cmap(200), legend.mar = 0.1,
             legend.width=0.8, legend.shrink=1,legend.cex = 0.6)
}

### GWR
for(l in 1:length(level.set))
{
  par(mfrow = c(2,2), oma = c(0, 0, 0, 3), mar = c(2, 1, 1, 2))
  zlim = c(-0.5,0.5)
  for(s in 1:n.season)
  {
    
    beta.long = data.frame(lon = ca_elevation$longitude, 
                           lat = ca_elevation$latitude,
                           b.gls = ca_tas_coef$ca_tas_nw$beta_est[,s],
                           b.gwr = ca_tas_coef$ca_tas_w$coef[2,s,]
    )
    se.long = data.frame(lon = ca_elevation$longitude, 
                         lat = ca_elevation$latitude,
                         b.gls = ca_tas_coef$ca_tas_nw$beta_se[,s],
                         b.gwr = ca_tas_coef$ca_tas_w$se.coef[2,s,]
    )
    
    beta.gls.wide = acast(beta.long, lon~lat, value.var='b.gls')
    beta.gwr.wide = acast(beta.long, lon~lat, value.var='b.gwr')
    se.gls.wide = acast(se.long, lon~lat, value.var='b.gls')
    se.gwr.wide = acast(se.long, lon~lat, value.var='b.gwr')
    mask = (!is.na(beta.gwr.wide))*1
    
    #Compute the residuals.
    R.gls= R.gwr = array(0,c(length(lon),length(lat),n))
    for(i in 1:n)
    {
      long.t = data.frame(lon = ca_elevation$longitude, 
                          lat = ca_elevation$latitude,
                          #R.gls = ks_tas_coef$ks_tas_nw$resids[,i,2],
                          R.gls = ca_tas_coef$ca_tas_nw$resids[,i,s],
                          R.gwr = ca_tas_coef$ca_tas_w$residuals[i,s,])
      R.gls[,,i] = acast(long.t, lon~lat, value.var='R.gls')
      R.gwr[,,i] = acast(long.t, lon~lat, value.var='R.gwr')
    }
    
    sigma.gls.hat = apply(R.gls,1:2,function(x) sd(x, na.rm = TRUE))
    R.gls.tilde = R.gls / rep(sigma.gls.hat,n)
    sigma.gwr.hat = apply(R.gwr,1:2,function(x) sd(x, na.rm = TRUE))
    R.gwr.tilde = R.gwr / rep(sigma.gwr.hat,n)
    
    sigma.gls.long = apply(ca_tas_coef$ca_tas_nw$resids[,,s],2,function(x) sd(x, na.rm = TRUE))
    R.gls.long = ca_tas_coef$ca_tas_nw$resids[,,s]/sigma.gls.long
    R.gls.sub = t(R.gls.long[select.ind,])
    
    sigma.gwr.long = apply(ca_tas_coef$ca_tas_w$residuals[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gwr.long = ca_tas_coef$ca_tas_w$residuals[,s,]/sigma.gwr.long
    R.gwr.sub = R.gwr.long[,select.ind]
    
    A = (mask==1)
    A[is.na(A)] = FALSE
    for(i in 1:n) R.gls.tilde[,,i] = R.gls.tilde[,,i] * A
    for(i in 1:n) R.gwr.tilde[,,i] = R.gwr.tilde[,,i] * A
    
    #Compute quantile (bootstrap).
    P_MC_gls = MC_gauss(R.gls.tilde,N)
    P_MC_gwr = MC_gauss(R.gwr.tilde,N)
    a_CoPE_gls = a_CoPE_gwr = 0
    while(P_MC_gls(a_CoPE_gls)>alpha_CoPE/3) a_CoPE_gls = a_CoPE_gls+0.01
    while(P_MC_gwr(a_CoPE_gwr)>alpha_CoPE/3) a_CoPE_gwr = a_CoPE_gwr+0.01
    
    #Normalized function.
    level= level.set[l]
    # norm_diff_gls = (beta.gls.wide - level) / se.gls.wide
    # norm_diff_gwr = (beta.gwr.wide - level) / se.gwr.wide
    
    n.lon = length(lon)
    n.lat = length(lat)
    
    norm_diff_gls = (beta.long$b.gls - level) / se.long$b.gls
    norm_diff_gwr = (beta.long$b.gwr - level) / se.long$b.gwr
    
    # image.map(lon,lat,beta.gls.wide,mask=mask,col = tim.colors(64),
    #           horizontal=FALSE,ylab='',xlab='',
    #           main= paste0("GLS ",season.var[s]," level=",level))
    
    gwr_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                norm_diff_gwr), 200, 200)$xyz.est
    mask_large = (!is.na(slope_interp$z))*1
    #zlim = c(0,0.6)
    slope_interp <- mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                   ca_tas_coef$ca_tas_w$coef[2,s,]), 200, 200)$xyz.est
    
    level = level.set[l]
    image.map(slope_interp$x,slope_interp$y,
              slope_interp$z,mask=mask_large,col = cmap(200),
              #horizontal=FALSE,
              legend = FALSE,
              ylab='',xlab='',
              zlim = zlim,
              main= paste0(season.var[s]))
    maps::map("county",add = TRUE)
    points(sub.loc$x, sub.loc$y,pch = 16, col = "black",cex=1.1)
    drawContour(slope_interp$x,slope_interp$y, bound = FALSE,
                slope_interp$z,c=level,col="purple",lty=1)
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=a_CoPE_gwr,col="darkred")
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=-a_CoPE_gwr,col="darkgreen")
    
  }
  par(mfrow=c(1, 1), mar=c(1, 1, 0.7, 2), new=FALSE)
  image.plot(zlim= zlim, legend.only=TRUE, col=cmap(200), legend.mar = 0.1,
             legend.width=0.8, legend.shrink=1,legend.cex = 0.6)
}

### future GWR
for(l in 1:length(level.set))
{
  par(mfrow = c(2,2), oma = c(0, 0, 0, 3), mar = c(2, 1, 1, 2))
  zlim = c(0,1.1)
  for(s in 1:n.season)
  {
    
    beta.long = data.frame(lon = ca_elevation$longitude, 
                           lat = ca_elevation$latitude,
                           b.gls = ca_tas_coef$ca_tas_nw$beta_est[,s],
                           b.gwr = ca_tas_coef$ca_tas_pre$coef[2,s,]
    )
    se.long = data.frame(lon = ca_elevation$longitude, 
                         lat = ca_elevation$latitude,
                         b.gls = ca_tas_coef$ca_tas_nw$beta_se[,s],
                         b.gwr = ca_tas_coef$ca_tas_pre$se.coef[2,s,]
    )
    
    beta.gls.wide = acast(beta.long, lon~lat, value.var='b.gls')
    beta.gwr.wide = acast(beta.long, lon~lat, value.var='b.gwr')
    se.gls.wide = acast(se.long, lon~lat, value.var='b.gls')
    se.gwr.wide = acast(se.long, lon~lat, value.var='b.gwr')
    mask = (!is.na(beta.gwr.wide))*1
    
    #Compute the residuals.
    R.gls= R.gwr = array(0,c(length(lon),length(lat),n))
    for(i in 1:n)
    {
      long.t = data.frame(lon = ca_elevation$longitude, 
                          lat = ca_elevation$latitude,
                          #R.gls = ks_tas_coef$ks_tas_nw$resids[,i,2],
                          R.gls = ca_tas_coef$ca_tas_nw$resids[,i,s],
                          R.gwr = ca_tas_coef$ca_tas_pre$residuals[i,s,])
      R.gls[,,i] = acast(long.t, lon~lat, value.var='R.gls')
      R.gwr[,,i] = acast(long.t, lon~lat, value.var='R.gwr')
    }
    
    sigma.gls.hat = apply(R.gls,1:2,function(x) sd(x, na.rm = TRUE))
    R.gls.tilde = R.gls / rep(sigma.gls.hat,n)
    sigma.gwr.hat = apply(R.gwr,1:2,function(x) sd(x, na.rm = TRUE))
    R.gwr.tilde = R.gwr / rep(sigma.gwr.hat,n)
    
    sigma.gls.long = apply(ca_tas_coef$ca_tas_nw$resids[,,s],2,function(x) sd(x, na.rm = TRUE))
    R.gls.long = ca_tas_coef$ca_tas_nw$resids[,,s]/sigma.gls.long
    R.gls.sub = t(R.gls.long[select.ind,])
    
    sigma.gwr.long = apply(ca_tas_coef$ca_tas_pre$residuals[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gwr.long = ca_tas_coef$ca_tas_pre$residuals[,s,]/sigma.gwr.long
    R.gwr.sub = R.gwr.long[,select.ind]
    
    A = (mask==1)
    A[is.na(A)] = FALSE
    for(i in 1:n) R.gls.tilde[,,i] = R.gls.tilde[,,i] * A
    for(i in 1:n) R.gwr.tilde[,,i] = R.gwr.tilde[,,i] * A
    
    #Compute quantile (bootstrap).
    P_MC_gls = MC_gauss(R.gls.tilde,N)
    P_MC_gwr = MC_gauss(R.gwr.tilde,N)
    a_CoPE_gls = a_CoPE_gwr = 0
    while(P_MC_gls(a_CoPE_gls)>alpha_CoPE/3) a_CoPE_gls = a_CoPE_gls+0.01
    while(P_MC_gwr(a_CoPE_gwr)>alpha_CoPE/3) a_CoPE_gwr = a_CoPE_gwr+0.01
    
    #Normalized function.
    level= level.set[l]
    # norm_diff_gls = (beta.gls.wide - level) / se.gls.wide
    # norm_diff_gwr = (beta.gwr.wide - level) / se.gwr.wide
    
    n.lon = length(lon)
    n.lat = length(lat)
    
    norm_diff_gls = (beta.long$b.gls - level) / se.long$b.gls
    norm_diff_gwr = (beta.long$b.gwr - level) / se.long$b.gwr
    
    # image.map(lon,lat,beta.gls.wide,mask=mask,col = tim.colors(64),
    #           horizontal=FALSE,ylab='',xlab='',
    #           main= paste0("GLS ",season.var[s]," level=",level))
    
    gwr_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                norm_diff_gwr), 200, 200)$xyz.est
    slope_interp <- mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                   ca_tas_coef$ca_tas_pre$coef[2,s,]), 200, 200)$xyz.est
    mask_large = (!is.na(slope_interp$z))*1
    #zlim = c(0,0.6)
    
    level = level.set[l]
    image.map(slope_interp$x,slope_interp$y,
              slope_interp$z,mask=mask_large,col = pos_cmap(200),
              #horizontal=FALSE,
              legend = FALSE,
              ylab='',xlab='',
              zlim = zlim,
              main= paste0(season.var[s]))
    maps::map("county",add = TRUE)
    points(sub.loc$x, sub.loc$y,pch = 16, col = "black",cex=1.1)
    drawContour(slope_interp$x,slope_interp$y, bound = FALSE,
                slope_interp$z,c=level,col="purple",lty=1)
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=a_CoPE_gwr,col="darkred")
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=-a_CoPE_gwr,col="darkgreen")
    
  }
  par(mfrow=c(1, 1), mar=c(1, 1, 0.7, 2), new=FALSE)
  image.plot(zlim= zlim, legend.only=TRUE, col=pos_cmap(200), legend.mar = 0.1,
             legend.width=0.8, legend.shrink=1,legend.cex = 0.6)
}

### GWR with new slopes
load("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/summary_results_20211101/ca_slope_v_new.rdata")
zlim = c(0,0.6)
level.set = c(0.15,0.2,0.25)
for(l in 1:length(level.set))
{
  par(mfrow = c(2,2), oma = c(0, 0, 0, 3), mar = c(2, 1, 1, 2))
  for(s in 1:n.season)
  {
    
    beta.long = data.frame(lon = ca_elevation$longitude, 
                           lat = ca_elevation$latitude,
                           b.gwr = slope_new$slope_hist$beta_w_tau[,s]
    )
    se.long = data.frame(lon = ca_elevation$longitude, 
                         lat = ca_elevation$latitude,
                         b.gwr = slope_new$slope_hist$sigma_w_tau[,s]
    )

    beta.gwr.wide = acast(beta.long, lon~lat, value.var='b.gwr')
    se.gwr.wide = acast(se.long, lon~lat, value.var='b.gwr')
    mask = (!is.na(beta.gwr.wide))*1
    
    #Compute the residuals.
    R.gwr = array(0,c(length(lon),length(lat),n))
    for(i in 1:n)
    {
      long.t = data.frame(lon = ca_elevation$longitude, 
                          lat = ca_elevation$latitude,
                          R.gwr = ca_tas_coef$ca_tas_w$residuals[i,s,])
      R.gwr[,,i] = acast(long.t, lon~lat, value.var='R.gwr')
    }
    
    sigma.gwr.hat = apply(R.gwr,1:2,function(x) sd(x, na.rm = TRUE))
    R.gwr.tilde = R.gwr / rep(sigma.gwr.hat,n)
    
    sigma.gwr.long = apply(ca_tas_coef$ca_tas_w$residuals[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gwr.long = ca_tas_coef$ca_tas_w$residuals[,s,]/sigma.gwr.long
    R.gwr.sub = R.gwr.long[,select.ind]
    
    A = (mask==1)
    A[is.na(A)] = FALSE
    for(i in 1:n) R.gwr.tilde[,,i] = R.gwr.tilde[,,i] * A
    
    #Compute quantile (bootstrap).
    P_MC_gwr = MC_gauss(R.gwr.tilde,N)
    a_CoPE_gwr = 0
    while(P_MC_gwr(a_CoPE_gwr)>alpha_CoPE/3) a_CoPE_gwr = a_CoPE_gwr+0.01
    
    #Normalized function.
    level= level.set[l]
    # norm_diff_gls = (beta.gls.wide - level) / se.gls.wide
    # norm_diff_gwr = (beta.gwr.wide - level) / se.gwr.wide
    
    n.lon = length(lon)
    n.lat = length(lat)
    
    norm_diff_gwr = (beta.long$b.gwr - level) / se.long$b.gwr
    
    # image.map(lon,lat,beta.gls.wide,mask=mask,col = tim.colors(64),
    #           horizontal=FALSE,ylab='',xlab='',
    #           main= paste0("GLS ",season.var[s]," level=",level))
    
    gwr_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                norm_diff_gwr), 200, 200)$xyz.est
    #zlim = c(0,0.6)
    mask_large = (!is.na(slope_interp$z))*1
    slope_interp <- mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                   slope_new$slope_hist$beta_w_tau[,s]), 200, 200)$xyz.est
    
    image.map(slope_interp$x,slope_interp$y,
              slope_interp$z,mask=mask_large,col = pos_cmap(200),
              #horizontal=FALSE,
              legend = FALSE,
              ylab='',xlab='',
              zlim = zlim,
              main= paste0(season.var[s]))
    maps::map("county",add = TRUE)
    points(sub.loc$x, sub.loc$y,pch = 16, col = "black",cex=1.1)
    drawContour(slope_interp$x,slope_interp$y, bound = FALSE,
                slope_interp$z,c=level,col="purple",lty=1)
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=a_CoPE_gwr,col="darkred")
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=-a_CoPE_gwr,col="darkgreen")
    
  }
  par(mfrow=c(1, 1), mar=c(1, 1, 0.7, 2), new=FALSE)
  image.plot(zlim= zlim, legend.only=TRUE, col=pos_cmap(200), legend.mar = 0.1,
             legend.width=0.8, legend.shrink=1,legend.cex = 0.6)
}

### Future GWR with new slopes
zlim=c(0,1)
level.set = c(0.5,0.55,0.6)
for(l in 1:length(level.set))
{
  par(mfrow = c(2,2), oma = c(0, 0, 0, 3), mar = c(2, 1, 1, 2))
  for(s in 1:n.season)
  {
    
    beta.long = data.frame(lon = ca_elevation$longitude, 
                           lat = ca_elevation$latitude,
                           b.gwr = slope_new$slope_rcp$beta_w_rcp_tau[,s]
    )
    se.long = data.frame(lon = ca_elevation$longitude, 
                         lat = ca_elevation$latitude,
                         b.gwr = slope_new$slope_rcp$sigma_w_rcp_tau[,s]
    )
    
    beta.gwr.wide = acast(beta.long, lon~lat, value.var='b.gwr')
    se.gwr.wide = acast(se.long, lon~lat, value.var='b.gwr')
    mask = (!is.na(beta.gwr.wide))*1
    
    #Compute the residuals.
    R.gwr = array(0,c(length(lon),length(lat),n))
    for(i in 1:n)
    {
      long.t = data.frame(lon = ca_elevation$longitude, 
                          lat = ca_elevation$latitude,
                          R.gwr = ca_tas_coef$ca_tas_pre$residuals[i,s,])
      R.gwr[,,i] = acast(long.t, lon~lat, value.var='R.gwr')
    }
    
    sigma.gwr.hat = apply(R.gwr,1:2,function(x) sd(x, na.rm = TRUE))
    R.gwr.tilde = R.gwr / rep(sigma.gwr.hat,n)
    
    sigma.gwr.long = apply(ca_tas_coef$ca_tas_pre$residuals[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gwr.long = ca_tas_coef$ca_tas_pre$residuals[,s,]/sigma.gwr.long
    R.gwr.sub = R.gwr.long[,select.ind]
    
    A = (mask==1)
    A[is.na(A)] = FALSE
    for(i in 1:n) R.gwr.tilde[,,i] = R.gwr.tilde[,,i] * A
    
    #Compute quantile (bootstrap).
    P_MC_gwr = MC_gauss(R.gwr.tilde,N)
    a_CoPE_gwr = 0
    while(P_MC_gwr(a_CoPE_gwr)>alpha_CoPE/3) a_CoPE_gwr = a_CoPE_gwr+0.01
    
    #Normalized function.
    level= level.set[l]
    # norm_diff_gls = (beta.gls.wide - level) / se.gls.wide
    # norm_diff_gwr = (beta.gwr.wide - level) / se.gwr.wide
    
    n.lon = length(lon)
    n.lat = length(lat)
    
    norm_diff_gwr = (beta.long$b.gwr - level) / se.long$b.gwr
    
    # image.map(lon,lat,beta.gls.wide,mask=mask,col = tim.colors(64),
    #           horizontal=FALSE,ylab='',xlab='',
    #           main= paste0("GLS ",season.var[s]," level=",level))
    
    gwr_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                norm_diff_gwr), 200, 200)$xyz.est
    #zlim = c(0,0.6)
    slope_interp <- mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                   slope_new$slope_rcp$beta_w_rcp_tau[,s]), 200, 200)$xyz.est
    mask_large = (!is.na(slope_interp$z))*1
    
    image.map(slope_interp$x,slope_interp$y,
              slope_interp$z,mask=mask_large,col = pos_cmap(200),
              #horizontal=FALSE,
              legend = FALSE,
              ylab='',xlab='',
              zlim = zlim,
              main= paste0(season.var[s]))
    maps::map("county",add = TRUE)
    points(sub.loc$x, sub.loc$y,pch = 16, col = "black",cex=1.1)
    drawContour(slope_interp$x,slope_interp$y, bound = FALSE,
                slope_interp$z,c=level,col="purple",lty=1)
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=a_CoPE_gwr,col="darkred")
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=-a_CoPE_gwr,col="darkgreen")
    
  }
  par(mfrow=c(1, 1), mar=c(1, 1, 0.7, 2), new=FALSE)
  image.plot(zlim= zlim, legend.only=TRUE, col=pos_cmap(200), legend.mar = 0.1,
             legend.width=0.8, legend.shrink=1,legend.cex = 0.6)
}

#### Precipitation
## 6.5 Cope sets

### GR
level.set = c(-0.05,0,0.05)
for(l in 1:length(level.set))
{
  zlim = zlim = c(-0.2,0.2)
  par(mfrow = c(2,2), oma = c(0, 0, 0, 3), mar = c(2, 1, 1, 2))
  for(s in 1:n.season)
  {
    n = dim(gwgr.pre.ca$res.gr$dat_resid)[1]
    beta.long = data.frame(lon = ca_elevation$longitude, 
                           lat = ca_elevation$latitude,
                           b.gls = gwgr.pre.ca$res.gr$gwgr.slope[,s],
                           b.gwr = gwgr.pre.ca$res.gwgr$gwgr.slope[,s]
    )
    se.long = data.frame(lon = ca_elevation$longitude, 
                         lat = ca_elevation$latitude,
                         b.gls = gwgr.pre.ca$res.gr$gwgr.slope.se[,s],
                         b.gwr = gwgr.pre.ca$res.gwgr$gwgr.slope.se[,s]
    )
    
    beta.gls.wide = acast(beta.long, lon~lat, value.var='b.gls')
    beta.gwr.wide = acast(beta.long, lon~lat, value.var='b.gwr')
    se.gls.wide = acast(se.long, lon~lat, value.var='b.gls')
    se.gwr.wide = acast(se.long, lon~lat, value.var='b.gwr')
    mask = (!is.na(beta.gwr.wide))*1
    
    #Compute the residuals.
    R.gls = R.gwr = array(0,c(length(lon),length(lat),n))
    for(i in 1:n)
    {
      long.t = data.frame(lon = ca_elevation$longitude, 
                          lat = ca_elevation$latitude,
                          #R.deviance = resid_mat[,s,i],
                          #R.gls = ks_tas_coef$ks_tas_nw$resids[,i,2],
                          R.gls = gwgr.pre.ca$res.gr$dat_resid[i,s,],
                          R.gwr = gwgr.pre.ca$res.gwgr$dat_resid[i,s,])
      #R.deviance[,,i] = acast(long.t, lon~lat, value.var='R.deviance')
      R.gls[,,i] = acast(long.t, lon~lat, value.var='R.gls')
      R.gwr[,,i] = acast(long.t, lon~lat, value.var='R.gwr')
    }
    
    sigma.gls.hat = apply(R.gls,1:2,function(x) sd(x, na.rm = TRUE))
    R.gls.tilde = R.gls / rep(sigma.gls.hat,n)
    sigma.gwr.hat = apply(R.gwr,1:2,function(x) sd(x, na.rm = TRUE))
    R.gwr.tilde = R.gwr / rep(sigma.gwr.hat,n)
    R.gls.sub = gwgr.pre.ca$res.gr$dat_resid[,s,select.ind]
    sigma.gls.long = apply(gwgr.pre.ca$res.gr$dat_resid[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gls.long = gwgr.pre.ca$res.gr$dat_resid[,s,]/sigma.gls.long
    R.gls.sub = R.gls.long[,select.ind]
    
    sigma.gwr.long = apply(gwgr.pre.ca$res.gwgr$dat_resid[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gwr.long = gwgr.pre.ca$res.gwgr$dat_resid[,s,]/sigma.gwr.long
    R.gwr.sub = R.gwr.long[,select.ind]
    A = (mask==1)
    A[is.na(A)] = FALSE
    for(i in 1:n) R.gls.tilde[,,i] = R.gls.tilde[,,i] * A
    for(i in 1:n) R.gwr.tilde[,,i] = R.gwr.tilde[,,i] * A
    
    #Compute quantile (bootstrap).
    P_MC_gls = MC_gauss(R.gls.tilde,N)
    P_MC_gwr = MC_gauss(R.gwr.tilde,N)
    a_CoPE_gls = a_CoPE_gwr = 0
    while(P_MC_gls(a_CoPE_gls)>alpha_CoPE) a_CoPE_gls = a_CoPE_gls+0.01
    while(P_MC_gwr(a_CoPE_gwr)>alpha_CoPE) a_CoPE_gwr = a_CoPE_gwr+0.01
    
    #Normalized function.
    level= level.set[l]
    # norm_diff_gls = (beta.gls.wide - level) / se.gls.wide
    # norm_diff_gwr = (beta.gwr.wide - level) / se.gwr.wide
    
    n.lon = length(lon)
    n.lat = length(lat)
    norm_diff_gls = (beta.long$b.gls - level) / se.long$b.gls
    norm_diff_gwr = (beta.long$b.gwr - level) / se.long$b.gwr
    
    # image.map(lon,lat,beta.gls.wide,mask=mask,col = tim.colors(64),
    #           horizontal=FALSE,ylab='',xlab='',
    #           main= paste0("GLS ",season.var[s]," level=",level))
    slope_interp <- mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                   gwgr.pre.ca$res.gr$gwgr.slope[,s]), 200, 200)$xyz.est
    gls_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                norm_diff_gls), 200, 200)$xyz.est
    gwr_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                norm_diff_gwr), 200, 200)$xyz.est
    mask_large = (!is.na(slope_interp$z))*1
    #zlim = range(-abs(slope_interp$z),abs(slope_interp$z),na.rm = TRUE)
    image.map(slope_interp$x,slope_interp$y,
              slope_interp$z,mask=mask_large,col = rev.cmap(200),
              #horizontal=FALSE,
              legend = FALSE,
              ylab='',xlab='',
              zlim = zlim,
              main= paste0(season.var[s]))
    maps::map("county",add = TRUE)
    points(sub.loc$x, sub.loc$y,pch = 16, col = "black",cex=1.1)
    drawContour(slope_interp$x,slope_interp$y,bound = FALSE,
                slope_interp$z,c=level,col="purple",lty=1)
    drawContour(gls_interp$x,gls_interp$y,gls_interp$z,c=a_CoPE_gls,col="darkgreen")
    drawContour(gls_interp$x,gls_interp$y,gls_interp$z,c=-a_CoPE_gls,col="darkred")
  }
  par(mfrow=c(1, 1), mar=c(1, 1, 0.7, 2), new=FALSE)
  image.plot(zlim= zlim, legend.only=TRUE, col=rev.cmap(200), legend.mar = 0.1,
             legend.width=0.8, legend.shrink=1,legend.cex = 0.6)
}


## GWGR

#par(mfrow = c(2,2), oma=c(1, 1, 1,1),  mar = c(2.5, 2.5, 2.5, 2.5))
for(l in 1:length(level.set))
{
  zlim = zlim = c(-0.2,0.2)
  par(mfrow = c(2,2), oma = c(0, 0, 0, 3), mar = c(2, 1, 1, 2))
  for(s in 1:n.season)
  {
    n = dim(gwgr.pre.ca$res.gr$dat_resid)[1]
    beta.long = data.frame(lon = ca_elevation$longitude, 
                           lat = ca_elevation$latitude,
                           b.gls = gwgr.pre.ca$res.gr$gwgr.slope[,s],
                           b.gwr = gwgr.pre.ca$res.gwgr$gwgr.slope[,s]
    )
    se.long = data.frame(lon = ca_elevation$longitude, 
                         lat = ca_elevation$latitude,
                         b.gls = gwgr.pre.ca$res.gr$gwgr.slope.se[,s],
                         b.gwr = gwgr.pre.ca$res.gwgr$gwgr.slope.se[,s]
    )
    
    beta.gls.wide = acast(beta.long, lon~lat, value.var='b.gls')
    beta.gwr.wide = acast(beta.long, lon~lat, value.var='b.gwr')
    se.gls.wide = acast(se.long, lon~lat, value.var='b.gls')
    se.gwr.wide = acast(se.long, lon~lat, value.var='b.gwr')
    mask = (!is.na(beta.gwr.wide))*1
    
    #Compute the residuals.
    R.gls = R.gwr = array(0,c(length(lon),length(lat),n))
    for(i in 1:n)
    {
      long.t = data.frame(lon = ca_elevation$longitude, 
                          lat = ca_elevation$latitude,
                          #R.deviance = resid_mat[,s,i],
                          #R.gls = ks_tas_coef$ks_tas_nw$resids[,i,2],
                          R.gls = gwgr.pre.ca$res.gr$dat_resid[i,s,],
                          R.gwr = gwgr.pre.ca$res.gwgr$dat_resid[i,s,])
      #R.deviance[,,i] = acast(long.t, lon~lat, value.var='R.deviance')
      R.gls[,,i] = acast(long.t, lon~lat, value.var='R.gls')
      R.gwr[,,i] = acast(long.t, lon~lat, value.var='R.gwr')
    }
    
    sigma.gls.hat = apply(R.gls,1:2,function(x) sd(x, na.rm = TRUE))
    R.gls.tilde = R.gls / rep(sigma.gls.hat,n)
    sigma.gwr.hat = apply(R.gwr,1:2,function(x) sd(x, na.rm = TRUE))
    R.gwr.tilde = R.gwr / rep(sigma.gwr.hat,n)
    R.gls.sub = gwgr.pre.ca$res.gr$dat_resid[,s,select.ind]
    sigma.gls.long = apply(gwgr.pre.ca$res.gr$dat_resid[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gls.long = gwgr.pre.ca$res.gr$dat_resid[,s,]/sigma.gls.long
    R.gls.sub = R.gls.long[,select.ind]
    
    sigma.gwr.long = apply(gwgr.pre.ca$res.gwgr$dat_resid[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gwr.long = gwgr.pre.ca$res.gwgr$dat_resid[,s,]/sigma.gwr.long
    R.gwr.sub = R.gwr.long[,select.ind]
    A = (mask==1)
    A[is.na(A)] = FALSE
    #for(i in 1:n) R.gls.tilde[,,i] = R.gls.tilde[,,i] * A
    for(i in 1:n) R.gwr.tilde[,,i] = R.gwr.tilde[,,i] * A
    
    #Compute quantile (bootstrap).
    P_MC_gls = MC_gauss(R.gls.tilde,N)
    P_MC_gwr = MC_gauss(R.gwr.tilde,N)
    a_CoPE_gls = a_CoPE_gwr = 0
    #while(P_MC_gls(a_CoPE_gls)>alpha_CoPE) a_CoPE_gls = a_CoPE_gls+0.01
    while(P_MC_gwr(a_CoPE_gwr)>alpha_CoPE/3) a_CoPE_gwr = a_CoPE_gwr+0.01
    
    #Normalized function.
    level= level.set[l]
    # norm_diff_gls = (beta.gls.wide - level) / se.gls.wide
    # norm_diff_gwr = (beta.gwr.wide - level) / se.gwr.wide
    
    n.lon = length(lon)
    n.lat = length(lat)
    norm_diff_gls = (beta.long$b.gls - level) / se.long$b.gls
    norm_diff_gwr = (beta.long$b.gwr - level) / se.long$b.gwr
    
    # image.map(lon,lat,beta.gls.wide,mask=mask,col = tim.colors(64),
    #           horizontal=FALSE,ylab='',xlab='',
    #           main= paste0("GLS ",season.var[s]," level=",level))
    gwr_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                norm_diff_gwr), 200, 200)$xyz.est
    slope_interp <- mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                   gwgr.pre.ca$res.gwgr$gwgr.slope[,s]), 200, 200)$xyz.est
    mask_large = (!is.na(slope_interp$z))*1
    image.map(slope_interp$x,slope_interp$y,
              slope_interp$z,mask=mask_large,col = rev.cmap(200),
              #horizontal=FALSE,
              ylab='',xlab='',
              zlim = zlim,
              legend = FALSE,
              main= paste0(season.var[s]))
    maps::map("county",add = TRUE)
    drawContour(slope_interp$x,slope_interp$y,bound=FALSE,
                slope_interp$z,c=level,col="purple",lty=1)
    points(sub.loc$x, sub.loc$y,pch = 16, col = "black",cex=1.1)
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=a_CoPE_gwr,col="darkgreen")
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=-a_CoPE_gwr,col="darkred")
  }
  par(mfrow=c(1, 1), mar=c(1, 1, 0.7, 2), new=FALSE)
  image.plot(zlim= zlim, legend.only=TRUE, col=rev.cmap(200), legend.mar = 0.1,
             legend.width=0.8, legend.shrink=1,legend.cex = 0.6)
}

## future GWGR
#par(mfrow = c(2,2), oma=c(1, 1, 1,1),  mar = c(2.5, 2.5, 2.5, 2.5))
for(l in 1:length(level.set))
{
  zlim = zlim = c(-0.22,0.22)
  par(mfrow = c(2,2), oma = c(0, 0, 0, 3), mar = c(2, 1, 1, 2))
  for(s in 1:n.season)
  {
    n = dim(gwgr.pre.ca$res.gr$dat_resid)[1]
    beta.long = data.frame(lon = ca_elevation$longitude, 
                           lat = ca_elevation$latitude,
                           b.gls = gwgr.pre.ca$res.gr$gwgr.slope[,s],
                           b.gwr = gwgr.pre.ca$future.gwgr$gwgr.slope[,s]
    )
    se.long = data.frame(lon = ca_elevation$longitude, 
                         lat = ca_elevation$latitude,
                         b.gls = gwgr.pre.ca$res.gr$gwgr.slope.se[,s],
                         b.gwr = gwgr.pre.ca$future.gwgr$gwgr.slope.se[,s]
    )
    
    beta.gls.wide = acast(beta.long, lon~lat, value.var='b.gls')
    beta.gwr.wide = acast(beta.long, lon~lat, value.var='b.gwr')
    se.gls.wide = acast(se.long, lon~lat, value.var='b.gls')
    se.gwr.wide = acast(se.long, lon~lat, value.var='b.gwr')
    mask = (!is.na(beta.gwr.wide))*1
    
    #Compute the residuals.
    R.gls = R.gwr = array(0,c(length(lon),length(lat),n))
    for(i in 1:n)
    {
      long.t = data.frame(lon = ca_elevation$longitude, 
                          lat = ca_elevation$latitude,
                          #R.deviance = resid_mat[,s,i],
                          #R.gls = ks_tas_coef$ks_tas_nw$resids[,i,2],
                          R.gls = gwgr.pre.ca$res.gr$dat_resid[i,s,],
                          R.gwr = gwgr.pre.ca$future.gwgr$dat_resid[i,s,])
      #R.deviance[,,i] = acast(long.t, lon~lat, value.var='R.deviance')
      R.gls[,,i] = acast(long.t, lon~lat, value.var='R.gls')
      R.gwr[,,i] = acast(long.t, lon~lat, value.var='R.gwr')
    }
    
    sigma.gls.hat = apply(R.gls,1:2,function(x) sd(x, na.rm = TRUE))
    R.gls.tilde = R.gls / rep(sigma.gls.hat,n)
    sigma.gwr.hat = apply(R.gwr,1:2,function(x) sd(x, na.rm = TRUE))
    R.gwr.tilde = R.gwr / rep(sigma.gwr.hat,n)
    R.gls.sub = gwgr.pre.ca$res.gr$dat_resid[,s,select.ind]
    sigma.gls.long = apply(gwgr.pre.ca$res.gr$dat_resid[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gls.long = gwgr.pre.ca$res.gr$dat_resid[,s,]/sigma.gls.long
    R.gls.sub = R.gls.long[,select.ind]
    
    sigma.gwr.long = apply(gwgr.pre.ca$future.gwgr$dat_resid[,s,],2,function(x) sd(x, na.rm = TRUE))
    R.gwr.long = gwgr.pre.ca$future.gwgr$dat_resid[,s,]/sigma.gwr.long
    R.gwr.sub = R.gwr.long[,select.ind]
    A = (mask==1)
    A[is.na(A)] = FALSE
    #for(i in 1:n) R.gls.tilde[,,i] = R.gls.tilde[,,i] * A
    for(i in 1:n) R.gwr.tilde[,,i] = R.gwr.tilde[,,i] * A
    
    #Compute quantile (bootstrap).
    P_MC_gls = MC_gauss(R.gls.tilde,N)
    P_MC_gwr = MC_gauss(R.gwr.tilde,N)
    a_CoPE_gls = a_CoPE_gwr = 0
    #while(P_MC_gls(a_CoPE_gls)>alpha_CoPE) a_CoPE_gls = a_CoPE_gls+0.01
    while(P_MC_gwr(a_CoPE_gwr)>alpha_CoPE/3) a_CoPE_gwr = a_CoPE_gwr+0.01
    
    #Normalized function.
    level= level.set[l]
    # norm_diff_gls = (beta.gls.wide - level) / se.gls.wide
    # norm_diff_gwr = (beta.gwr.wide - level) / se.gwr.wide
    
    n.lon = length(lon)
    n.lat = length(lat)
    norm_diff_gls = (beta.long$b.gls - level) / se.long$b.gls
    norm_diff_gwr = (beta.long$b.gwr - level) / se.long$b.gwr
    
    # image.map(lon,lat,beta.gls.wide,mask=mask,col = tim.colors(64),
    #           horizontal=FALSE,ylab='',xlab='',
    #           main= paste0("GLS ",season.var[s]," level=",level))
    gwr_interp = mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                norm_diff_gwr), 200, 200)$xyz.est
    slope_interp <- mba.surf(cbind(ca_elevation[,c("longitude","latitude")],
                                   gwgr.pre.ca$future.gwgr$gwgr.slope[,s]), 200, 200)$xyz.est
    mask_large = (!is.na(slope_interp$z))*1
    image.map(slope_interp$x,slope_interp$y,
              slope_interp$z,mask=mask_large,col = rev.cmap(200),
              #horizontal=FALSE,
              ylab='',xlab='',
              zlim = zlim,
              legend = FALSE,
              main= paste0(season.var[s]))
    maps::map("county",add = TRUE)
    drawContour(slope_interp$x,slope_interp$y,bound=FALSE,
                slope_interp$z,c=level,col="purple",lty=1)
    points(sub.loc$x, sub.loc$y,pch = 16, col = "black",cex=1.1)
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=a_CoPE_gwr,col="darkgreen")
    drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=-a_CoPE_gwr,col="darkred")
  }
  par(mfrow=c(1, 1), mar=c(1, 1, 0.7, 2), new=FALSE)
  image.plot(zlim= zlim, legend.only=TRUE, col=rev.cmap(200), legend.mar = 0.1,
             legend.width=0.8, legend.shrink=1,legend.cex = 0.6)
}
dev.off()