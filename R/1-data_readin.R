rm(list=ls())
library(reshape2)
library(ncdf4)
library(sf)

data_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/raw/pre_tas/"
setwd(data_path)
all_files = list.files()
nco = ncdf4::nc_open(all_files[1])
lon = ncdf4::ncvar_get(nco, varid = "lon")
lat = ncdf4::ncvar_get(nco, varid = "lat")
ncdf4::nc_close(nco)

# determine which coords are in each state polygon
map_path = "/Users/wenyilin/Documents/R/NA-CORDEX/map/"
within_ca = readRDS(file = paste0(map_path,"within_ca.rds"))
within_co = readRDS(file = paste0(map_path,"within_co.rds"))
within_ks = readRDS(file = paste0(map_path,"within_ks.rds"))
ca_poly = readRDS(file = paste0(map_path,"ca_poly.rds"))
co_poly = readRDS(file = paste0(map_path,"co_poly.rds"))
ks_poly = readRDS(file = paste0(map_path,"ks_poly.rds"))

# create sf coords of lon/lat
dtf = expand.grid(lon, lat)
names(dtf) = c("lon", "lat")
coords = sf::st_as_sf(dtf, coords = c("lon", "lat"),
                      crs = sf::st_crs(ks_poly))
# save data in .rds
save_path = "/Users/wenyilin/Documents/R/NA-CORDEX/data/rds/pre_tas/"
for (i in seq_along(all_files)) {
  # open connection and load relevant data
  nco = ncdf4::nc_open(all_files[i])
  message(paste(i, "/", length(all_files), ": ", all_files[i], sep = ""))
  lon = ncdf4::ncvar_get(nco, varid = "lon")
  lat = ncdf4::ncvar_get(nco, varid = "lat")
  coords = expand.grid(lon, lat)
  tas = ncdf4::ncvar_get(nco, varid = "pr")
  # close connection
  ncdf4::nc_close(nco)
  
  # extract relevant data for each state
  #tas = matrix(c(tas) - 273.15, ncol = dim(tas)[3]) # col: time
  tas = matrix(c(tas), ncol = dim(tas)[3]) 
  coords_ca = coords[within_ca,]
  tas_ca = tas[within_ca, ]
  coords_co = coords[within_co,]
  tas_co = tas[within_co, ]
  coords_ks = coords[within_ks,]
  tas_ks = tas[within_ks, ]
  
  # save data
  
  #saveRDS(coords_ca, file = paste0(save_path,"coords_ca_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  saveRDS(tas_ca, file = paste0(save_path,"pr_ca_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  #saveRDS(coords_co, file = paste0("coords_co_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  saveRDS(tas_co, file = paste0(save_path,"pr_co_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  #saveRDS(coords_ks, file = paste0("coords_ks_", all_files[i], ".rds", sep = ""), compress = "bzip2")
  saveRDS(tas_ks, file = paste0(save_path,"pr_ks_", all_files[i], ".rds", sep = ""), compress = "bzip2")
}
