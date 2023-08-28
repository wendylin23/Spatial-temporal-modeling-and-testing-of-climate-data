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
library(Vizumap)
library(nlme)
library(autoimage)
library(FRK)
library(pwr)
library(MTS)
library(tseries)
library(RColorBrewer)
library(dglm)
library(quantreg)
code_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/R/"
data_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Data" 
res_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Results"
setwd(code_path)
source("varx_fixed.R")
source("util.R")

#####################################
############## Load data ############
#####################################
#tas_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/tas-rec-rcp85-mon-44i/"
#pr_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/pr-rec-rcp85-mon-44i/"
map_path = paste0(data_path,'/map/')
within_ca = readRDS(file = paste0(map_path,"within_ca.rds"))
load(paste0(map_path,"within_rec_ca.rdata"))
load(paste0(map_path,"ca_elevation_rec.rdata"))
load(paste0(map_path,"ca_geodata.rdata"))
## color map
cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
rev.cmap <- colorRampPalette(brewer.pal(11, "RdBu"))
neg_cmap = colorRampPalette(rev(brewer.pal(11, "RdBu"))[1:6]) ## all blue
pos_cmap = colorRampPalette(rev((brewer.pal(11, "RdBu"))[1:6])) ## all red
## color map for data
pr.cmap <- colorRampPalette(brewer.pal(11, "Spectral")[6:11])
tas.cmap <- colorRampPalette(rev(brewer.pal(11, "Spectral"))[1:6])
ele.cmap <- colorRampPalette(terrain.colors(10, alpha = 1))


#################################################
##### Load and prepare temperature data #########
#################################################
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

## exploratory plots
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

#################################################
##### Load and prepare precipitation data #######
#################################################

pr_hist_ca = readRDS(file = paste0(data_path,'/pre_tas/pr_ca_prec.hist.CanESM2.CanRCM4.mon.NAM-44i.raw.nc.rds'))
pr_rcp_ca = readRDS(file = paste0(data_path,'/pre_tas/pr_ca_prec.rcp85.CanESM2.CanRCM4.mon.NAM-44i.raw.nc.rds'))

## seasonal precipitation data
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

ca_data = list(gls.tas.summary = gls.tas.summary,
               gls.tas.rcp = gls.tas.rcp,
               gls.pr.summary = gls.pr.summary,
               gls.pr.rcp = gls.pr.rcp,
               ca_elevation = ca_elevation)
# save(ca_data, file = "yourdatapath/ca_data.rdata") ## Save CA temperature and precipitation data
