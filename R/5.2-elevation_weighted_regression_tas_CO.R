rm(list=ls())
library(forecast)
library(reshape2)
library(ggplot2)
library(ggmap)
#library(maps)
#library(mapdata)
library(mgcv)
library(lme4)
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
#library(spgwr)

## get aspect ratio from lon/lat
map_aspect = function(x, y) {
  x.center <- sum(range(x)) / 2
  y.center <- sum(range(y)) / 2
  x.dist <- ggplot2:::dist_central_angle(x.center + c(-0.5, 0.5), rep(y.center, 2))
  y.dist <- ggplot2:::dist_central_angle(rep(x.center, 2), y.center + c(-0.5, 0.5))
  y.dist / x.dist
}

tas_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/tas-rcp85-mon-44i/"
pr_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/pr-rcp85-mon-44i/"
map_path = "/Users/wenyilin/Documents/R/NA-CORDEX/map/"
code_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/"
res_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/"
#within_co = readRDS(file = paste0(map_path,"within_co.rds"))
#load(paste0(map_path,"within_co.rdata"))
load(paste0(map_path,"co_elevation.rdata"))
## Get elevation
#Sys.setenv(id = "wendylin23")
#user = Sys.getenv("id")
#co_elevation = elevation(longitude = coords_co[,1],
                         # latitude = coords_co[,2],elevation_model = "srtm1",
                         # username = user, verbose = FALSE)
tas_hist_co_all = list.files(path = tas_path,pattern = "^tas_co_tas.hist")
tas_rcp_co_all = list.files(path = tas_path,pattern = "^tas_co_tas.rcp85")
### example analysis for CA
tas_hist_co = readRDS(file = paste0(tas_path,tas_hist_co_all[1]))
tas_rcp_co = readRDS(file = paste0(tas_path,tas_rcp_co_all[1]))
co_elevation$ele_std = (co_elevation$elevation_geonames - min(co_elevation$elevation_geonames))/(max(co_elevation$elevation_geonames)-min(co_elevation$elevation_geonames))
co_elevation$ele_norm = (co_elevation$elevation_geonames - mean(co_elevation$elevation_geonames))/1000

## seasonal tas data
time_dat = data.frame(year = rep(1950:2005,each=12),
                      month = rep(1:12,56),
                      ts.point = 1:(12*56))
gls.dat = data.frame(lat = co_elevation[,1], lon = co_elevation[,2],
                     elevation = co_elevation$ele_norm,
                     tas = tas_hist_co,
                     loc = as.factor(1:length(co_elevation[,1])))
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
gls.dat = data.frame(lat = co_elevation[,1], lon = co_elevation[,2],
                     elevation = co_elevation$ele_norm,
                     tas = tas_rcp_co,
                     loc = as.factor(1:length(co_elevation[,1])))
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

## basic definitions
season.var = c("winter","spring","summer","fall" )
loc.names = c("Denver","Aspen","Grand Junction")
select.loc = data.frame(dv = c(39.75,-104.75), ap = c(39.25,-106.75),
                        gj = c(39.25,-108.75))
select.ind = c(79,61,57)
select.loc.col = c("red","green","blue")
sub.loc = list(x = c(-104.75,-106.75,-108.75), 
               y = c(39.75,39.25,39.25),
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
## color map
cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
rev.cmap <- colorRampPalette(brewer.pal(11, "RdBu"))
neg_cmap = colorRampPalette(rev(brewer.pal(11, "RdBu"))[1:6]) ## all blue
pos_cmap = colorRampPalette(rev((brewer.pal(11, "RdBu"))[1:6])) ## all red

## spatial seasonal trend
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
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
          zlim = c(-12,32),
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
autoimage(co_elevation$longitude,co_elevation$latitude,
          co_elevation$elevation_geonames,
          interp.args = list(no.X = 200, no.Y = 200),
          common.legend = FALSE,
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = c("elevation (m) in CO"),
          # outer.title = "",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          #size = c(2,3), lratio = 0.4,
          #col.args = list(cmap1(100),cmap2(100))``
          #col=cmap1(100),
          #zlim = list(c(10,35),c(0,3500)),
          legend = "vertical")

## weighted regression on elevation
source(paste0(code_path,"varx_fixed.R"))

beta_est = beta_se = matrix(0, nrow = n.loc, ncol = n.season)
b0_est = matrix(0, nrow = n.loc, ncol = n.season)
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
  b0_est[i,] = m_ar0_cor2$coef[1,]
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
autoimage(co_elevation$longitude,co_elevation$latitude,
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
autoimage(co_elevation$longitude,co_elevation$latitude,
          z_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.25,
          col=cmap(200),
          zlim = c(-3.5,3.5),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          season_coef,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=pos_cmap(200),
          zlim = c(0,0.7),
          legend = "vertical")

z_score = season_coef/season_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          z_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.25,
          col=pos_cmap(200),
          zlim = c(0,5),
          legend = "vertical")

co_tas_nw = list(b0_est = b0_est,
                 beta_est = beta_est,
                 beta_se = beta_se,
                 season_coef = season_coef,
                 season_se = season_se,
                 resids = resids)

## test spatial correlation
rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) 
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}
n.lon = length(unique(gls.tas.summary$lon))
n.lat = length(unique(gls.tas.summary$lat))
simgrid = expand.grid(1:n.lon, 1:n.lat)
distance = as.matrix(dist(simgrid))
r.distance = as.matrix(dist(simgrid)) * 50 ##0.44 degree = 50km

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

par(mfrow=c(2,2))
for(t in 1:n.season)
{
  tmp.cor = as.vector(cor.list[[t]])
  vec.dist = as.vector(distance)
  tmp.dat = data.frame(tmp.cor = tmp.cor, vec.dist = vec.dist)
  plot(vec.dist, tmp.cor, pch=".", cex=0.8, main = season.var[t])
  abline(lm(tmp.cor ~ vec.dist, data = tmp.dat), col="red")
  tmp.quad = lm(tmp.cor ~ vec.dist + I(vec.dist^2), data = tmp.dat)
  dist = seq(0, max(vec.dist), by=0.1)
  pred.cor = predict(tmp.quad,newdata = data.frame(vec.dist = dist))
  lines(dist, pred.cor, col="blue")
  #legend("topright", legend = c("linear","quadratic"),col = c("red","blue"), lty=1)
}

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
  #pred.cor = predict(fit_lm2,newdata = data.frame(vec.dist = dist.seq))
  pred.cor = predict(fit_lm21,newdata = data.frame(dist2 = dist.seq^2))
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

## special case: wii = 1 and wij=0 if i\neq j
source(paste0(code_path,"gwr.R"))
m1 = lm(summer ~ time.pts, data = gls.tas.summary)
m2.1 = lm(winter ~ time.pts + elevation, data = gls.tas.summary)
m2.2 = lm(spring ~ time.pts + elevation + winter, data = gls.tas.summary)
m2.21 = lmer(spring ~ time.pts + elevation + winter + (1|loc), data = gls.tas.summary)
m2.3 = lm(summer ~ time.pts + elevation + spring, data = gls.tas.summary)
m2.31 = lmer(summer ~ time.pts + elevation + spring + (1|loc), data = gls.tas.summary)
m2.4 = lm(fall ~ time.pts + elevation+ summer, data = gls.tas.summary)

w0 = diag(1, n.loc)
tas_exp = array(0, dim = c(n.time-1,4,n.loc))
tas_rcp = array(0, dim = c(n.rcp.time-1,4,n.loc))
x_mult = array(0, dim = c(n.time-1,7,n.loc))
x_rcp = array(0, dim = c(n.rcp.time-1,7,n.loc))
fixed_beta1 = rbind(matrix(1,nrow = 2,ncol = 4),diag(1,4))
for(i in 1:n.loc)
{
  data.i = gls.tas.summary[gls.tas.summary$loc==i,]
  data.rcp.i = gls.tas.rcp[gls.tas.rcp$loc==i,]
  tas_exp[,,i] = as.matrix(data.i[-1,season.var])
  tas_rcp[,,i] = as.matrix(data.rcp.i[-1,season.var])
  x_mult[,,i] =  cbind(rep(1,n.time-1),
                       time.pts[-1],
                       rep(data.i$elevation[1],n.time-1), #with elevation
                       data.i$fall[1:(n.time-1)],data.i$winter[2:n.time],
                       data.i$spring[2:n.time],data.i$summer[2:n.time])
  x_rcp[,,i] =  cbind(rep(1,n.rcp.time-1),
                      time.rcp.pts[-1],
                      rep(data.rcp.i$elevation[1],n.rcp.time-1), #with elevation
                      data.rcp.i$fall[1:(n.rcp.time-1)],data.rcp.i$winter[2:n.rcp.time],
                      data.rcp.i$spring[2:n.rcp.time],data.rcp.i$summer[2:n.rcp.time])
}
fit_m1 = varw_fixed(tas_exp,p=0,xt = x_mult[,-3,],w = w0,fixed = fixed_beta1) ## no elevation
fit_m2 = varw_fixed(tas_exp,p=0,xt = x_mult[,-1,],w = w0,fixed = fixed_beta1) ## with elevation
w_beta_est = t(fit_m2$coef[1,,])
w_beta_se = t(fit_m2$se.coef[1,,])

w_ele_est = t(fit_m2$coef[2,,])
w_ele_se = t(fit_m2$se.coef[2,,])

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_beta_est,
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

z_score = w_beta_est/w_beta_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          z_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.25,
          col=cmap(200),
          zlim = c(-3.5,3.5),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_ele_est,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (elevation) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=cmap(200),
          zlim = c(-0.03,0.03),
          legend = "vertical")

ele_z_score = w_ele_est/w_ele_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          ele_z_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of estimated slopes (elevation) of four seasons",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.25,
          col=cmap(200),
          zlim = c(-100,100),
          legend = "vertical")

## bisquare kernel
source(paste0(code_path,"gwr.R"))
d_bis = seq(1,6,by = 0.5)

### Model selection via AIC
fixed_beta2 = rbind(matrix(1,nrow = 3,ncol = 4),diag(1,4))
fixed_beta3 = rbind(matrix(1,nrow = 3,ncol = 4),diag(0,4))
aic.bisquare.w1 = aic.bisquare.w2 = rss.w2 = matrix(0, nrow = 4, ncol = length(d_bis))
for(j in 1:length(d_bis))
{
  w1 = gwr.bisquare(distance^2,d_bis[j])
  w2 = gwr.lm(distance,d_bis[j])
  fit_t1 = varw_fixed(tas_exp,p=0,xt = x_mult,w = w1,fixed = fixed_beta2)
  fit_t2 = varw_fixed(tas_exp,p=0,xt = x_mult,w = w2,fixed = fixed_beta2)
  aic.bisquare.w1[,j] = fit_t1$aic
  aic.bisquare.w2[,j] = fit_t2$aic
  rss.w2[,j] = fit_t2$rss
}## d=1.6

d_gauss = seq(1,10,by = 0.5)
aic.gauss.w3 = aic.gauss.w4 = matrix(0, nrow = 4, ncol = length(d_gauss))
for(j in 1:length(d_gauss))
{
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
  for(j in 1:length(d_bis))
  {
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
  }## d=1.6
}

source(paste0(code_path,"gwr.R"))
fixed_beta2 = rbind(matrix(1,nrow = 3,ncol = 4),diag(1,4))
fixed_beta3 = rbind(matrix(1,nrow = 3,ncol = 4),diag(0,4))
w1 = gwr.bisquare(distance^2,4)
w2 = gwr.lm(distance,2.5)
w3 = gwr.exp(distance, 2)
w4 = gwr.gauss(distance^2,2)

fit_m3_2 = varw_fixed(tas_exp,p=0,xt = x_mult,w = w2,fixed = fixed_beta2) ## selected weight matrix
### Fit the model with future data
pred_m3_2 = varw_fixed(tas_rcp,p=0,xt = x_rcp,w = w2,fixed = fixed_beta2)

# co_tas_coef = list(co_tas_nw = co_tas_nw,
#                    co_tas_w = fit_m3_2,
#                    co_tas_pre = pred_m3_2)
# save(co_tas_coef,file = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/summary_results_20211101/co_tas_coef.rdata")

w_aic = rbind(fit_m3_1$aic,fit_m3_2$aic,fit_m3_3$aic,fit_m3_4$aic)
w_b0_est = list(t(fit_m3_1$coef[1,,]),t(fit_m3_2$coef[1,,]),
                t(fit_m3_3$coef[1,,]),t(fit_m3_4$coef[1,,]))
w_slope_est = list(t(fit_m3_1$coef[2,,]),t(fit_m3_2$coef[2,,]),
                   t(fit_m3_3$coef[2,,]),t(fit_m3_4$coef[2,,]))
w_h_est = list(t(fit_m3_1$coef[3,,]),t(fit_m3_2$coef[3,,]),
               t(fit_m3_3$coef[3,,]),t(fit_m3_4$coef[3,,]))
w_sea_est = list(cbind(fit_m3_1$coef[4,1,],fit_m3_1$coef[5,2,],
                       fit_m3_1$coef[6,3,],fit_m3_1$coef[7,4,]),
                 cbind(fit_m3_2$coef[4,1,],fit_m3_2$coef[5,2,],
                       fit_m3_2$coef[6,3,],fit_m3_2$coef[7,4,]),
                 cbind(fit_m3_3$coef[4,1,],fit_m3_3$coef[5,2,],
                       fit_m3_3$coef[6,3,],fit_m3_3$coef[7,4,]),
                 cbind(fit_m3_4$coef[4,1,],fit_m3_4$coef[5,2,],
                       fit_m3_4$coef[6,3,],fit_m3_4$coef[7,4,]))

w_b0_se = list(t(fit_m3_1$se.coef[1,,]),t(fit_m3_2$se.coef[1,,]),
               t(fit_m3_3$se.coef[1,,]),t(fit_m3_4$se.coef[1,,]))
w_slope_se = list(t(fit_m3_1$se.coef[2,,]),t(fit_m3_2$se.coef[2,,]),
                  t(fit_m3_3$se.coef[2,,]),t(fit_m3_4$se.coef[2,,]))
w_h_se = list(t(fit_m3_1$se.coef[3,,]),t(fit_m3_2$se.coef[3,,]),
              t(fit_m3_3$se.coef[3,,]),t(fit_m3_4$se.coef[3,,]))
w_sea_se = list(cbind(fit_m3_1$se.coef[4,1,],fit_m3_1$se.coef[5,2,],
                      fit_m3_1$se.coef[6,3,],fit_m3_1$se.coef[7,4,]),
                cbind(fit_m3_2$se.coef[4,1,],fit_m3_2$se.coef[5,2,],
                      fit_m3_2$se.coef[6,3,],fit_m3_2$se.coef[7,4,]),
                cbind(fit_m3_3$se.coef[4,1,],fit_m3_3$se.coef[5,2,],
                      fit_m3_3$se.coef[6,3,],fit_m3_3$se.coef[7,4,]),
                cbind(fit_m3_4$se.coef[4,1,],fit_m3_4$se.coef[5,2,],
                      fit_m3_4$se.coef[6,3,],fit_m3_4$se.coef[7,4,]))

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_b0_est[[3]],
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated intercept (C) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=cmap(200),
          #zlim = c(-10,32),
          zlim = c(-ceiling(max(abs(w_b0_est[[1]]))),ceiling(max(abs(w_b0_est[[1]])))),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_slope_est[[3]],
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slope (C/decade) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=cmap(200),
          #zlim = c(-0.32,0.32),
          zlim = c(-0.3,0.3),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_h_est[[3]],
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated elevation slope (C/km) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=neg_cmap(100),
          #zlim = c(-0.32,0.32),
          zlim = c(-3.5,0),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_sea_est[[3]],
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated seasonal correaltion (C/C) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=pos_cmap(100),
          #zlim = c(-0.32,0.32),
          zlim = c(0,1),
          legend = "vertical")

## z-scores
z_score = list(slope_m1 = w_slope_est[[1]]/w_slope_se[[1]],
               slope_m2 = w_slope_est[[2]]/w_slope_se[[2]],
               slope_m3 = w_slope_est[[3]]/w_slope_se[[3]],
               slope_m4 = w_slope_est[[4]]/w_slope_se[[4]])
z_score_nw = beta_est/beta_se

co_res = list(cor.rsquare = cor.rsquare,
              aic.list = aic.list,
              beta_nw = beta_est,
              z_score_nw = z_score_nw, w_b0_est = w_b0_est, 
              w_slope_est=w_slope_est,w_h_est=w_h_est,
              w_sea_est = w_sea_est,z_score_w = z_score)
save(co_res,file = paste0(res_path,"co_res.rdata"))

par(mfrow=c(2,2),mar=c(2,2,3,1))
col1 = rgb(1,0,0,1/4)
col2 = rgb(0,1,0,1/4)
col3 = rgb(0,0,1,1/4)
for(i in 1:n.season)
{
  tmp.z = cbind(z_score$slope_m2[,i],z_score$slope_m3[,i],z_score_nw[,i])
  hist(tmp.z[,1],xlab = "Z-score",main = paste0(season.var[i]," slope z scores"), breaks = 30, col=col1,xlim = range(tmp.z))
  #hist(z_score$slope_m1[,2],col=rgb(0,1,1,1/4), add = T, breaks =30)
  hist(tmp.z[,2],col=col2, add = T, breaks =30)
  #hist(z_score$slope_m4[,2],col=rgb(1,1,0,1/4), add = T, breaks =30)
  hist(tmp.z[,3],col=col3, add = T, breaks =30)
  if(i==n.season)
    legend("topleft",legend = c("lm_weighted","exp_weighted","non-weighted"),
         fill = c(col1,col2,col3),
         bty = 'n',cex = 0.75,
         border = NA)
}

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          z_score$slope_m1,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z scores of time slope of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=cmap(200),
          #zlim = c(-0.32,0.32),
          zlim = c(-6,6),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          z_score$h_z,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z scores of elevation slope of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=cmap(200),
          #zlim = c(-0.32,0.32),
          zlim = c(-10,10),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(ca_elevation$longitude,ca_elevation$latitude,
          z_score$sea_z,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z scores of seasonal correlation of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=cmap(200),
          #zlim = c(-0.32,0.32),
          zlim = c(0,18),
          legend = "vertical")

## comparing estimated slopes
cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_beta_est,
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
          zlim = c(-0.32,0.32),
          legend = "vertical")
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_beta_est2,
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
          zlim = c(-0.32,0.32),
          legend = "vertical")

## comparing z_score
z_score1 = w_beta_est/w_beta_se
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          z_score1,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.25,
          col=cmap(200),
          zlim = c(-4.5,4.5),
          legend = "vertical")

z_score2 = w_beta_est2/w_beta_se2
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          z_score2,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of estimated slopes of four seasons",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.25,
          col=cmap(200),
          zlim = c(-4.5,4.5),
          legend = "vertical")

par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_ele_est2,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Estimated slopes (elevation) of four seasons",
          mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          legend.axis.args = list(cex.axis=1),
          size = c(2, 2), lratio = 0.25,
          col=cmap(200),
          zlim = c(-0.03,0.03),
          legend = "vertical")

z_score = w_ele_est2/w_ele_se2
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          z_score,
          # interp.args = list(no.X = 200, no.Y = 200),
          map = "county",ylab = "Latitude",xlab = "Longitude",
          main = season.var,
          outer.title = "Z_scores of estimated slopes (elevation) of four seasons",
          mtext.args = list(cex = 1),
          legend.axis.args = list(cex.axis=1),
          size = c(2,2), lratio = 0.25,
          col=cmap(200),
          zlim = c(-80,80),
          legend = "vertical")

## Gaussian kernel
w4 = gwr.Gauss(distance^2,d[3])
fit_m4_1 = varw_fixed(tas_exp,p=0,xt = x_mult[,-3,],w = w4,fixed = fixed_beta1)
fit_m4_2 = varw_fixed(tas_exp,p=0,xt = x_mult[,-1,],w = w4,fixed = fixed_beta1)
fit_m4_3 = varw_fixed(tas_exp,p=0,xt = x_mult[,-3,],w = w4,fixed = fixed_beta2)
fit_m4_4 = varw_fixed(tas_exp,p=0,xt = x_mult[,-1,],w = w4,fixed = fixed_beta2)
w_beta_est4 = t(fit_m4_4$coef[1,,])
w_beta_se4 = t(fit_m4$se.coef[1,,])
w_ele_est4 = t(fit_m4$coef[2,,])
w_ele_se4 = t(fit_m4$se.coef[2,,])

## comparing estimated slopes
cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_beta_est4,
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
          zlim = c(-0.32,0.32),
          legend = "vertical")
par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
autoimage(co_elevation$longitude,co_elevation$latitude,
          w_beta_est2,
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
          zlim = c(-0.32,0.32),
          legend = "vertical")

##### pacakge checking
gwr.bl = gwr.sel(summer ~ time.pts, gls.tas.summary,
                 coords = cbind(gls.tas.summary$lon,gls.tas.summary$lat),
                 adapt = TRUE, longlat = TRUE)
gwr.bl2 = gwr(summer ~ time.pts, gls.tas.summary,verbose = TRUE,
              coords = as.matrix(gls.tas.summary$elevation),
              adapt = TRUE, longlat = TRUE)
gwr.model = gwr(Leaf.ec~Year,data=chinensis2,
                bandwidth=750000,fit.points=us_grid2)

