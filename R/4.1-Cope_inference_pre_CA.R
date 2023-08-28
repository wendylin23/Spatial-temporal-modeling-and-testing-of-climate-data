rm(list=ls())
library(grDevices)
library(cope)
library(reshape2)
library(autoimage)
library(fields)
library(maps)
library(RColorBrewer)
library(car)
code_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/R/"
data_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Data" 
res_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Results"
setwd(code_path)
source(paste0(code_path,"gwr.R"))
source(paste0(code_path,"varx_fixed.R"))
source(paste0(code_path,"util.R"))

## color map for parameters
cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
rev.cmap <- colorRampPalette(brewer.pal(11, "RdBu"))
neg_cmap = colorRampPalette(rev(brewer.pal(11, "RdBu"))[1:6]) ## all blue
pos_cmap = colorRampPalette(rev((brewer.pal(11, "RdBu"))[1:6])) ## all red

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

## Hist - pre
Simul_Cope(map = ca_elevation, par = gwgr.pre.ca$res.gwgr$gwgr.slope,
           par.se =  gwgr.pre.ca$res.gwgr$gwgr.slope.se,
           resids = gwgr.pre.ca$res.gwgr$dat_resid,
           level.set = c(-0.05,0,0.05), zlim = c(-0.2,0.2),
           col.map = rev.cmap(200),type = "prep",
           sub.loc=sub.loc, select.ind = select.ind)

## RCP - pre
Simul_Cope(map = ca_elevation, par = gwgr.pre.ca$future.gwgr$gwgr.slope,
           par.se =  gwgr.pre.ca$future.gwgr$gwgr.slope.se,
           resids = gwgr.pre.ca$future.gwgr$dat_resid,
           level.set = c(-0.05,0,0.05), zlim = c(-0.25,0.25),
           col.map = rev.cmap(200),type = "prep",
           sub.loc=sub.loc, select.ind = select.ind)
