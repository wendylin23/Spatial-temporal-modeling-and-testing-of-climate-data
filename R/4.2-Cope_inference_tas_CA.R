###########################################################################################
## This code file provides the analysis pipeline for deriving COPE set for temperature data
###########################################################################################

rm(list=ls())
library(grDevices)
library(cope)
library(reshape2)
library(autoimage)
library(fields)
library(maps)
library(RColorBrewer)
library(car)
code_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/"
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

## RCP - temp
Simul_Cope(map = ca_elevation, par = slope_new$slope_rcp$beta_w_rcp_tau,
           par.se =  slope_new$slope_rcp$sigma_w_rcp_tau,
           resids = ca_tas_coef$ca_tas_pre$residuals,
           level.set = c(0.5,0.55,0.6), zlim = c(0.4,1),
           col.map = pos_cmap(200),
           sub.loc=sub.loc, select.ind = select.ind)