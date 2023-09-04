#######################################################################################
## This code file provide the source code to read and transform downloaded data
## Raw Data can be downloaded from https://www.earthsystemgrid.org/search/cordexsearch.html
## Check the box in each column to select subset specified in this work
## CLick 'Download Options for Selection'
## Click 'Download Curl Script' and run the script from command line to install .nc data
## Save the downloaded data to local directory and using the below code to extract and 
## tranform data at state level.
#######################################################################################

rm(list=ls())
library(reshape2)
library(ncdf4)
library(sf)

# Specify dataset path
data_path = "/Users/wenyilin/Documents/GitHub/Spatial-temporal-modeling-and-testing-of-climate-data/Data" 
setwd(data_path)
# all_files = list.files(path = paste0(data_path,'/map/'))
# nco = ncdf4::nc_open(path = paste0(data_path,'/map/'))
# lon = ncdf4::ncvar_get(nco, varid = "lon")
# lat = ncdf4::ncvar_get(nco, varid = "lat")
# ncdf4::nc_close(nco)

# determine which coords are in each state polygon
map_path = paste0(data_path,'/map/')
load(paste0(map_path,"within_rec_ca.rdata"))
within_ca = readRDS(file = paste0(map_path,"within_ca.rds"))
within_co = readRDS(file = paste0(map_path,"within_co.rds"))
within_ks = readRDS(file = paste0(map_path,"within_ks.rds"))
ca_poly = readRDS(file = paste0(map_path,"ca_poly.rds"))
co_poly = readRDS(file = paste0(map_path,"co_poly.rds"))
ks_poly = readRDS(file = paste0(map_path,"ks_poly.rds"))

# create sf coords of lon/lat
# dtf = expand.grid(lon, lat)
# names(dtf) = c("lon", "lat")
# coords = sf::st_as_sf(dtf, coords = c("lon", "lat"),
#                       crs = sf::st_crs(ks_poly))
# save data in .rds
all_files = list.files(path = paste0(data_path,'/climate_data/pr/')) # change the path to temperature data if working on it
setwd(paste0(data_path,'/climate_data/pr'))
save_path = paste0(data_path,'/pre_tas/')
for (i in seq_along(all_files)) {
  # open connection and load relevant data
  nco = ncdf4::nc_open(all_files[i])
  message(paste(i, "/", length(all_files), ": ", all_files[i], sep = ""))
  lon = ncdf4::ncvar_get(nco, varid = "lon")
  lat = ncdf4::ncvar_get(nco, varid = "lat")
  coords = expand.grid(lon, lat)
  tmp = ncdf4::ncvar_get(nco, varid = "prec") # varid = 'temp' if referring to temperature data
  # close connection
  ncdf4::nc_close(nco)
  
  # extract relevant data for each state
  #tas = matrix(c(tas) - 273.15, ncol = dim(tas)[3]) # col: time
  tmp = matrix(c(tmp), ncol = dim(tmp)[3]) 
  coords_ca = coords[within_rec_ca,]
  tmp_ca = tmp[within_rec_ca, ]
  coords_co = coords[within_co,]
  tmp_co = tmp[within_co, ]
  coords_ks = coords[within_ks,]
  tmp_ks = tmp[within_ks, ]
  
  # save data
  ### Start with 'tas' if reffering to temperature data
  #saveRDS(coords_ca, file = paste0(save_path,"coords_ca_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  saveRDS(tmp_ca, file = paste0(save_path,"pr_ca_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  #saveRDS(coords_co, file = paste0("coords_co_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  saveRDS(tmp_co, file = paste0(save_path,"pr_co_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  #saveRDS(coords_ks, file = paste0("coords_ks_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  saveRDS(tmp_ks, file = paste0(save_path,"pr_ks_", all_files[i], ".rds", sep = ""), compress = "bzip2")
}
