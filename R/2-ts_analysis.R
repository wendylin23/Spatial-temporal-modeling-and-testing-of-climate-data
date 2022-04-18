rm(list = ls())
library(reshape2)
library(mgcv)
library(TSA)
library(dplyr)
library(fda)

data_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/tas-rcp85-mon-44i/"
map_path = "/Users/wenyilin/Documents/R/NA-CORDEX/map/"
#within_ca = readRDS(file = paste0(map_path,"within_ca.rds"))
#within_co = readRDS(file = paste0(map_path,"within_co.rds"))
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

### time series analysis
time_dat = data.frame(year = rep(1950:2005,each=12),
                      month = rep(1:12,56))

tas_hist_ts = ts(tas_hist_ks[1,],start=1950,frequency=12)
ts.plot(tas_hist_ts, ylab="Temperature")  

tas_exp = tas_hist_ts
time.pts = c(1:length(tas_exp))
time.pts = c(time.pts - min(time.pts))/max(time.pts)
mav.fit = ksmooth(time.pts, tas_exp, kernel = "box") ## fit with moving average
temp.fit.mav = ts(mav.fit$y,start=1950,frequency=12)  
lm.fit = lm(tas_exp~time.pts) ## fit with linear regression
temp.fit.lm = ts(fitted(lm.fit),start=1950,frequency=12)
loc.fit = loess(tas_exp~time.pts)  ## fit with local polynomial 
temp.fit.loc = ts(fitted(loc.fit),start=1950,frequency=12)
gam.fit = gam(tas_exp~s(time.pts))  ## fit with splines
temp.fit.gam = ts(fitted(gam.fit),start=1950,frequency=12) 

###plot the trend
ts.plot(tas_exp,ylab="Temperature")
lines(temp.fit.mav,lwd=2,col="purple")
lines(temp.fit.lm,lwd=2,col="green")
lines(temp.fit.loc,lwd=2,col="brown")
lines(temp.fit.gam,lwd=2,col="red")
abline(temp.fit.mav[1],0,lwd=2,col="blue")

### compare all results
all.val = c(temp.fit.mav,temp.fit.lm,temp.fit.gam,temp.fit.loc)
ylim= c(min(all.val),max(all.val))
ts.plot(temp.fit.lm,lwd=2,col="green",ylim=ylim,ylab="Temperature")
lines(temp.fit.mav,lwd=2,col="purple")
lines(temp.fit.gam,lwd=2,col="red")
lines(temp.fit.loc,lwd=2,col="brown")
legend("topleft",legend=c("MAV","LM","GAM","LOESS"),
       lty = 1, col=c("purple","green","red","brown"))

## estimate with Fourier basis
month = season(tas_exp)
har2 = harmonic(tas_exp,2)
temp.fit.harl = lm(tas_exp ~ har2)
summary(temp.fit.harl)
temp.fit.harl.gam = gam(tas_exp~s(time.pts)+har2)

dif.fit.harl = ts((tas_exp-fitted(temp.fit.harl)),start=1950,frequency=12)
dif.fit.harl.gam = ts((tas_exp-fitted(temp.fit.harl.gam)),start=1950,frequency=12)
ts.plot(dif.fit.harl,ylab="Residual Process")
lines(dif.fit.harl.gam,col="blue")

## Seasonal analysis
tas_hist_ks_long = cbind(time_dat,t(tas_hist_ks))
tas_hist_ks_mon = data.frame(tas_hist_ks_long[,-1] %>% group_by(month) %>% summarise_all("mean")
)
tas_hist_ks_mon_long = melt(tas_hist_ks_mon, id.vars=c("month"))
plot(tas_hist_ks_mon_long$month,tas_hist_ks_mon_long$value,
     pch = 1, cex = 0.5, col='lightblue',
     xlab="month", ylab="mean temperature")
fit.gam <- gam(value~s(month, k = 10, bs = "cr"),
               data = tas_hist_ks_mon_long, method="REML")
yhat <- predict(fit.gam,newdata = data.frame(month=1:12))
plot(tas_hist_ks_mon_long$month,tas_hist_ks_mon_long$value,
     pch = 1, cex = 0.5, col='lightblue',
     xlab="month", ylab="mean temperature")
lines(1:12, yhat,lwd=2, col="black")

## Spatial-Temporal Analysis with mgcv
coords_ks$ind = 1:nrow(coords_ks)
coords_wide = acast(coords_ks, Var1~Var2, value.var="ind")
select_spot = diag(coords_wide[nrow(coords_wide):1, ])
tas_sub_ks = cbind(time_dat,t(tas_hist_ks)[,select_spot])
tas_hist_ks_all_long = melt(tas_sub_ks, id.vars=c("year","month"))
fbasis=create.fourier.basis(rangeval = c(1,672), nbasis=5, period=12)
bvals = eval.basis(tas_hist_ks_all_long$month,fbasis)
fit.fourier = lm(tas_hist_ks_all_long$value ~ bvals[,-1]+as.factor(tas_hist_ks_all_long$variable))
tas_hist_ks_all_long$fitted = as.vector(fit.fourier$fitted.values)
tas_ks_exp_fitted = data.frame(tas_hist_ks_all_long[,-1] %>% group_by(month,variable) %>% summarize_all("mean"))
plot(tas_ks_exp_fitted$month, tas_ks_exp_fitted$value, pch = 1, cex = 0.5, col=tas_ks_exp_fitted$variable,
     xlab="month", ylab="temperature")
for(i in 1:6)
{
  lines(data=tas_ks_exp_fitted[tas_ks_exp_fitted$variable==i,],fitted ~ month, lwd=2, col=i)   # fitted
}

plot(1:nrow(tas_hist_ks_all_long), tas_hist_ks_all_long$value, pch = 1, cex = 0.5, col=tas_hist_ks_all_long$variable,
     xlab="month", ylab="temperature")
for(i in 1:6)
{
  lines(data=tas_ks_exp_fitted[tas_ks_exp_fitted$variable==i,],fitted ~ 1:nrow(tas_ks_exp_fitted), lwd=2, col=i)   # fitted
}



