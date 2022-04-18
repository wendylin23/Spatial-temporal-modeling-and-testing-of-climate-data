#####################################
########## Load packages ############
#####################################
library(forecast)
library(reshape2)
library(ggplot2)
library(mgcv)
library(TSA)
library(dplyr)
library(fda)
library(sp)
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
code_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/"
setwd(code_path)
source("varx_fixed.R")
source("util.R")

#####################################
############## Load data ############
#####################################
tas_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/tas-rec-rcp85-mon-44i/"
pr_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/pr-rec-rcp85-mon-44i/"
map_path = "/Users/wenyilin/Documents/R/NA-CORDEX/map/"
code_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/"
res_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/"
load(paste0(res_path,"summary_results_20211101/ks_data.rdata"))

#####################################
############ Load data ##############
#####################################

time_dat = data.frame(year = rep(1950:2005,each=12),
                     month = rep(1:12,56),
                     ts.point = 1:(12*56))  ## historical measures
time_rcp = data.frame(year = rep(2006:2100,each=12),
                      month = rep(1:12,95),
                      ts.point = 1:(12*95))
## temperature
gls.tas.summary = ks_data$gls.tas.summary
gls.tas.rcp = ks_data$gls.tas.rcp

## precipitation
gls.pr.summary = ks_data$gls.pr.summary
gls.pr.rcp = ks_data$gls.pr.rcp

## location and elevation
ks_elevation = ks_data$ks_elevation

#####################################
###########  plotting ###############
#####################################

## basic definitions
season.var = c("winter","spring","summer","fall" )
loc.names = c("Kansas City","Wichita","Oakley")
select.loc = data.frame(kc = c(39.25,-94.75), wi = c(37.75,-97.25),
                        oa = c(39.25,-100.75))
select.ind = c(75,25,63)
select.loc.col = c("red","green","blue")
sub.loc = list(x = c(-94.97,-97.25,-100.75), 
               y = c(39.25,37.75,39.25),
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
gls.tas.all = rbind(gls.tas.summary,gls.tas.summary)

## color map
cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
rev.cmap <- colorRampPalette(brewer.pal(11, "RdBu"))
neg_cmap = colorRampPalette(rev(brewer.pal(11, "RdBu"))[1:6]) ## all blue
pos_cmap = colorRampPalette(rev((brewer.pal(11, "RdBu"))[1:6])) ## all red

#png("figures/tas_4season_map.png",width = 250, height = 100, res = 100,units = "in")
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ks_elevation$longitude,ks_elevation$latitude,
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
          size = c(2, 2), lratio = 0.25,
          legend = "vertical")

## temporal seasonal trend
### historical trends
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
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,cex = 1.1,
       legend=season.var,lty = 1, col=1:length(season.var))
mtext("Average seasonal temperature (C) at three typical locations", outer = TRUE, cex = 1.2)

#####################################
######### TS model selection ########
#####################################

# Temperature
## 1. normal test with shapiro test
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
  plot_list(d.resids, main = season.var[j])
}

normal.res = table(norm_test < 0.05)
names(normal.res) = c("normal","non-normal")
normal.res

## 2. ACF and PACF
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

boxplot(lm_pacf, ylim = c(-acf_sig-0.1,acf_sig+0.1))
abline(h = acf_sig,col="red",lty=2)
abline(h = -acf_sig,col="red",lty=2)

x = box_pval[order(box_pval)]
plot(ecdf(box_pval))
abline(0,1)
lines(x, x +2*sqrt(x*(1-x)/90))
lines(x, x -2*sqrt(x*(1-x)/90))

## 3. cross-correlation between seasons
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
boxplot(multi.pval)
abline(h=0.05,col="red")

## 4. Test with AIC
aic_all = matrix(0,nrow = n.loc, ncol = 6)
for(i in 1:n.loc)
{
  tas_exp = gls.tas.summary[gls.tas.summary$loc==i,season.var]
  m_ar1 = varx_fixed(tas_exp,p=1,xt = matrix(time.pts, nrow = length(time.pts)))
  fixed_beta = rbind(rep(1,4),
                     diag(1,4),
                     rep(1,4))
  m_ar1_1 = varx_fixed(tas_exp[-1,],p=1,
                       xt = matrix(time.pts, nrow = length(time.pts)),
                       fixed = fixed_beta)
  m_ar0 = varx_fixed(tas_exp,p=0,xt = matrix(time.pts, nrow = length(time.pts)))
  fixed_beta = rbind(matrix(1,nrow = 6,ncol = 4),
                     diag(1,4))
  m_ar1_cor = varx_fixed(tas_exp[-1,],p=1,xt = x_multi,fixed = fixed_beta)
  fixed_beta = rbind(rep(1,4),
                     diag(1,4),
                     rep(1,4),
                     diag(1,4))
  m_ar1_nocor = varx_fixed(tas_exp[-1,],p=1,xt = x_multi,fixed = fixed_beta)
  fixed_beta = rbind(matrix(1,nrow = 2,ncol = 4),diag(1,4))
  m_ar0_cor = varx_fixed(tas_exp[-1,],p=0,xt = x_multi,fixed = fixed_beta)
  aic_all[i,] = round(c(m_ar1$aic,m_ar1_1$aic,m_ar0$aic,m_ar1_cor$aic,m_ar1_nocor$aic,m_ar0_cor$aic),4)
}
colnames(aic_all) = c("AR1 + t","nocorAR1 + t","t", "AR1 + t + season_lag","nocorAR1 + t +season_lag","t+season_lag")

# Precipitation
## 1. ACF and PACF
lm_acf = lm_pacf = matrix(NA, nrow = n.loc, ncol = n.season)
dw_p = dw_val = matrix(NA, nrow = n.loc, ncol = n.season)
resid_mat = array(0, dim = c(n.loc, n.season, n.time))
acf_sig = qnorm(1 - 0.05 / 2) / sqrt(n.time)

for(i in 1:n.loc)
{
  for(j in 1:n.season)
  {
    tas_exp = gls.pr.summary[gls.pr.summary$loc==i,season.var[j]]
    glm.fit = glm(tas_exp ~ time.pts, family=Gamma(link = "log"))
    resid_mat[i,j,] = residuals(glm.fit, type = c("deviance"))
    dw_res = durbinWatsonTest.glm(glm.fit)
    dw_p[i,j] = dw_res$p
    dw_val[i,j] = dw_res$dw
    acf.res = acf(residuals(glm.fit),plot = FALSE)
    lm_acf[i,j] = acf.res$acf[1,,]
    pacf.res = pacf(residuals(glm.fit),plot = FALSE)
    lm_pacf[i,j] = pacf.res$acf[1,,]
  }
}

dw_all = list(dw_p = dw_p, dw_val = dw_val)
#save(dw_all,file = paste0(res_path,"ks_pre_dw_p.rdata"))
boxplot(lm_pacf, ylim = c(-acf_sig-0.1,acf_sig+0.1))
abline(h = acf_sig,col="red",lty=2)
abline(h = -acf_sig,col="red",lty=2)

## 3. cross-correlation between seasons
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
boxplot(multi.pval)
abline(h=0.05,col="red")

#####################################
######### TS model fitting ##########
#####################################

### temperature
center_apply <- function(x) {
  apply(x, 2, function(y) y - mean(y))
}
beta_est = beta_se = matrix(0, nrow = n.loc, ncol = n.season)
season_coef = season_se = matrix(0, nrow = n.loc, ncol = n.season)
sresids = resids = array(NA, dim = c(n.loc, n.time-1, n.season))
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
  beta_se[i,] = m_ar0_cor2$se.beta[,1]
  season_coef[i,] = c(m_ar0_cor2$beta[1,2],m_ar0_cor2$beta[2,3],
                      m_ar0_cor2$beta[3,4],
                      m_ar0_cor2$beta[4,4])
  season_se[i,] = c(m_ar0_cor2$se.beta[1,2],m_ar0_cor2$se.beta[2,3],
                    m_ar0_cor2$se.beta[3,4],
                    m_ar0_cor2$se.beta[4,4])
  resids[i,,] = m_ar0_cor2$residuals
  #resids_c = center_apply(m_ar0_cor2$residuals)
  sresids[i,,] = t((chol(solve(m_ar0_cor2$Sigma)))%*% t(m_ar0_cor2$residuals))
}

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ks_elevation$longitude,ks_elevation$latitude,
          beta_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=cmap(200),
          zlim = c(-0.3,0.3),
          legend = "vertical")

z_score = beta_est/beta_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ks_elevation$longitude,ks_elevation$latitude,
          z_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.25,
          col=cmap(200),
          zlim = c(-2.5,2.5),
          legend = "vertical")

### precipitation
slope_est = slope_se =  p.value = gof.test = matrix(0, nrow = n.loc, ncol = n.season)
sig.level = 0.1
for(i in 1:n.loc)
{
  for(j in 1:n.season)
  {
    tas_exp = gls.pr.summary[gls.pr.summary$loc==i,season.var[j]]
    fit1 <- glm(tas_exp ~ time.pts, family=Gamma(link = "log"))
    deviance = fit1$deviance
    p.value[i,j] = pchisq(deviance, df = fit1$df.residual, lower.tail = F)
    gof.test[i,j] = ifelse(p.value[i,j] < sig.level,0,1)
    shape_est <- 1/summary(fit1)$dispersion
    scale_est <- 1/fit1$coefficients[1]/shape_est
    slope_est[i,j] <- fit1$coefficients[2]
    slope_se[i,j] <- summary(fit1)$coefficients[2,2]
  }
}

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ks_elevation$longitude,ks_elevation$latitude,
          slope_est,
          #interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (mm/month/10 years) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          col=rev.cmap(100),
          zlim = c(-0.1,0.1),
          size = c(2,2), lratio = 0.25,
          legend = "vertical")

z_score = slope_est/slope_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ks_elevation$longitude,ks_elevation$latitude,
          z_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z-scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=2),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=rev.cmap(100),
          zlim = c(-2.5,2.5),
          legend = "vertical")
