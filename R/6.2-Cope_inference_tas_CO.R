rm(list=ls())
library(grDevices)
library(cope)
library(reshape2)
library(autoimage)
library(fields)
library(maps)
library(RColorBrewer)
library(car)
source("~/Dropbox/UCSD/Thesis/3.Precipitation/Code/Max/functions.r")
## color map for parameters
cmap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
rev.cmap <- colorRampPalette(brewer.pal(11, "RdBu"))
neg_cmap = colorRampPalette(rev(brewer.pal(11, "RdBu"))[1:6]) ## all blue
pos_cmap = colorRampPalette(rev((brewer.pal(11, "RdBu"))[1:6])) ## all red


### Define functions
MBContour = function(x, y, R, cont, N = 1000){
  
  if(length(cont)==0) return(rep(-Inf,N))
  
  cont_x = unlist(sapply(cont, function(q) q$x))
  cont_y = unlist(sapply(cont, function(q) q$y))
  cont = cbind(cont_x,cont_y)
  
  interp_max = function(G){
    G = matrix(G,ncol=length(y))
    max(fields::interp.surface(list(x=x,y=y,z=G),cont), na.rm = TRUE)
  }
  
  n = dim(R)[3]
  g = matrix(rnorm(n*N),n,N)
  apply(abs(matrix(R,ncol=n) %*% g),2,interp_max) / sqrt(n-2)
}

get.rej2 = function(a_MB,mu_hat,R)
{
  a = a_MB
  Acz = mu_hat
  zna = R <= a & R >= -a
  which_na = which(zna == TRUE, arr.ind = TRUE)
  zmask = matrix(1, nrow = nrow(mu_hat), ncol = ncol(mu_hat))
  zmask[which_na] = NA
  z = Acz * zmask
  return(z)
}

plot_list = function(x, main = "", xlab = ""){
  maxy = sapply(x, function(j) j$y)
  plot(x[[1]], ylim = range(maxy), xlab = xlab, main = main, xlim = c(-3,3),
       col = "darkgray", lwd = 0.5)
  for (j in seq_along(x)) {
    lines(x[[j]], col = "darkgray", lwd = 0.5)
  }
}

MC_gauss = function(Y,N){
  n = dim(Y)[3]
  Y = matrix(Y,ncol=n)
  g = matrix(rnorm(n*N),n,N)
  maxima = apply(abs(Y %*% g), 2, function(x) max(x, na.rm = TRUE)) / sqrt(n)
  
  function(t) sum(maxima>=t) / length(maxima)
}

image.map = function(lon, lat, img, mask=NULL, xlab='longitude', ylab='latitude', ...) {
  if(!is.null(mask)) {
    img = img*mask
    xlim = lon[range(which(rowSums(mask, na.rm=TRUE)>0))]
    ylim = lat[range(which(colSums(mask, na.rm=TRUE)>0))]
  } else {
    xlim = range(lon)
    ylim = range(lat)
  }
  image.plot(lon, lat, img, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  #map('world', add=TRUE)
}

drawContour = function(x,y,z,c,col,lty=1){
  if(c>=0 & all(z>c))
  {
    #print("Full Set")
    text(max(x)-1, max(y)-1, "Full Upper Set", xpd=NA, col=col)
  }else if(c<=0 & all(z<c))
  {
    text(min(x)+1, min(y)+1, "Full Lower Set", xpd=NA, col=col)
  }else if(c>=0 & all(z<c))
  {
    text(max(x)-1, max(y)-1, "Empty Upper Set", xpd=NA, col=col)
  }else if(c<=0 & all(z>c))
  {
    #print("Empty Set")
    text(min(x)+1, min(y)+1, "Empty Lower Set", xpd=NA, col=col)
  }else
  {
    C<-contourLines(x,y,z,levels=c,nlevels=1)
    
    for(i in 1:length(C))
      lines(C[[i]]$x,C[[i]]$y,pch=20,col=col,lwd=3,lty=lty)
  }
}

qqnorm_new = function(y, pch = 16, xlab = "Standard Normal Quantiles",
                      add = FALSE, cex = 0.6, conf = 0.95,
                      ylab = "Sample Quantiles",
                      make.plot=TRUE,...)
{
  args = list(...)
  y = sort(na.omit(y))
  n <- length(y)
  P <- stats::ppoints(n)
  z <- qnorm(P)
  Q.x = quantile(y, c(.25,.75))
  Q.z = qnorm(c(.25,.75))
  b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
  a <- Q.x[1] - b*Q.z[1]
  
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/dnorm(z))*sqrt(P*(1 - P)/n)
  fit.value <- a + b*z
  upper <- fit.value + zz*SE
  lower <- fit.value - zz*SE
  if(make.plot & !add) {
    if(is.null(args$xlim)) plot(z, y, type="n", xlim = range(lower, z, upper, na.rm = TRUE), xlab = xlab, ylab = ylab, pch = pch,...)
    else plot(z, y, type="n", xlab = xlab, ylab = ylab, pch = pch, cex = cex, ...)
    points(z, y,pch=pch, cex=cex)
    abline(a, b, col = "blue", lwd=1)
    lines(z, upper, lty=2, lwd=0.8, col="blue")
    lines(z, lower, lty=2, lwd=0.8, col="blue")
  }else if(add)
  {
    points(z, y,pch=pch, cex=cex)
    abline(a, b, col = "blue", lwd=1)
    lines(z, upper, lty=2, lwd=0.8, col="blue")
    lines(z, lower, lty=2, lwd=0.8, col="blue")
  }
  out = data.frame(lower=lower, upper=upper, qnorm=z, data=y)
}

## Example for kansas
res_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/summary_results_20211101/"
code_path = "/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/"
load(paste0(res_path,"co_data.rdata"))
load(paste0(res_path,"co_pre_coef.rdata"))
load(paste0(res_path,"co_tas_coef.rdata"))
load(paste0(res_path,"co_tas_cv.rdata"))
load(paste0(res_path,"co_tas_aic_w2.rdata"))

season.var = c("winter","spring","summer","fall" )
n.season = length(season.var)
time_dat = data.frame(year = rep(1950:2005,each=12),
                      month = rep(1:12,56),
                      ts.point = 1:(12*56))  ## historical measures
time_rcp = data.frame(year = rep(2006:2100,each=12),
                      month = rep(1:12,95),
                      ts.point = 1:(12*95))
## temperature
gls.tas.summary = co_data$gls.tas.summary
gls.tas.rcp = co_data$gls.tas.rcp

## precipitation
gls.pr.summary = co_data$gls.pr.summary
gls.pr.rcp = co_data$gls.pr.rcp

## location and elevation
co_elevation = co_data$co_elevation

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
gls.tas.all = rbind(gls.tas.summary,gls.tas.summary)

gls.pr.summary$time.pts = time.pts
gls.pr.rcp$time.pts = time.rcp.pts
gls.pr.all = rbind(gls.pr.summary,gls.pr.rcp)

#### Cope set for Kansas
#Confidence level for CoPE sets.
alpha_CoPE = 0.1
#Nominal expected false area ratio for FARE sets.
alpha_FARE = 0.05
#Number of realizations to produce in the MC simulations.
N=1000
lon = unique(co_elevation$longitude)
lat = unique(co_elevation$latitude)
n = dim(co_tas_coef$co_tas_w$residuals)[1]

# Temperature data
for(s in 1:n.season)
{
  beta.long = data.frame(lon = co_elevation$longitude, 
                         lat = co_elevation$latitude,
                         b.gls = co_tas_coef$co_tas_nw$beta_est[,s],
                         b.gwr = co_tas_coef$co_tas_w$coef[2,s,]
  )
  se.long = data.frame(lon = co_elevation$longitude, 
                       lat = co_elevation$latitude,
                       b.gls = co_tas_coef$co_tas_nw$beta_se[,s],
                       b.gwr = co_tas_coef$co_tas_w$se.coef[2,s,]
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
    long.t = data.frame(lon = co_elevation$longitude, 
                        lat = co_elevation$latitude,
                        #R.gls = ks_tas_coef$ks_tas_nw$resids[,i,2],
                        R.gls = co_tas_coef$co_tas_nw$resids[,i,s],
                        R.gwr = co_tas_coef$co_tas_w$residuals[i,s,])
    R.gls[,,i] = acast(long.t, lon~lat, value.var='R.gls')
    R.gwr[,,i] = acast(long.t, lon~lat, value.var='R.gwr')
  }
  
  sigma.gls.hat = apply(R.gls,1:2,sd)
  R.gls.tilde = R.gls / rep(sigma.gls.hat,n)
  sigma.gwr.hat = apply(R.gwr,1:2,sd)
  R.gwr.tilde = R.gwr / rep(sigma.gwr.hat,n)
  
  sigma.gls.long = apply(gwgr.pre.co$res.gr$dat_resid[,s,],2,function(x) sd(x, na.rm = TRUE))
  R.gls.long = gwgr.pre.co$res.gr$dat_resid[,s,]/sigma.gls.long
  R.gls.sub = R.gls.long[,select.ind]
  
  sigma.gwr.long = apply(gwgr.pre.co$res.gwgr$dat_resid[,s,],2,function(x) sd(x, na.rm = TRUE))
  R.gwr.long = gwgr.pre.co$res.gwgr$dat_resid[,s,]/sigma.gwr.long
  R.gwr.sub = R.gwr.long[,select.ind]
  
  ### Assumptions for CoPE set
  ## A1. The estimated regression coefficients are continuous across space.
  par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
  autoimage(co_elevation$longitude,co_elevation$latitude,
            co_tas_coef$co_tas_nw$beta_est,
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
  
  par(oma = c(0, 0, 3, 0), mar = c(2, 1, 1, 1))
  autoimage(co_elevation$longitude,co_elevation$latitude,
            t(co_tas_coef$co_tas_w$coef[2,,]),
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
  
  ## A2. The estimated regression coefﬁcients are approximately Gaussian.
  norm_gls_test = norm_gwr_test = matrix(0, nrow = length(lon), ncol = length(lat))
  
  #plot(NULL, xlim=c(-5,5), ylim=c(-5,5), ylab="Observed", xlab="norm quantile")
  for(i in 1:length(lon))
  {
    for(j in 1:length(lat))
    {
      #points(R.gls.tilde[i,j,],pch=16,cex = 0.8, col=alpha("darkgrey",0.8))
      norm_gls_test[i,j] = shapiro.test(R.gls.tilde[i,j,])$p.value
      # d.resids = density(R.gls.tilde[i,j,])
      # if(i==1&j==1)
      #   plot(d.resids, ylim = c(0,0.8), xlab = "", main = "", xlim = c(-3,3),
      #        col = "darkgray", lwd = 0.5)
      # else
      #   lines(d.resids, col = "darkgray", lwd = 0.5)
      
      ## change to qqplot
    }
  }
  #lines(d.resids$x,dnorm(d.resids$x))
  
  qqnorm_new(R.gls.tilde[1,1,])
  
  #normal.res = table(norm_gls_test < 0.05)
  #normal.res
  
  for(i in 1:length(lon))
  {
    for(j in 1:length(lat))
    {
      norm_gwr_test[i,j] = shapiro.test(R.gwr.tilde[i,j,])$p.value
    }
  }
  
  #normal.res = table(norm_gwr_test < 0.05)
  #normal.res
  
  ##A3. The transformed errors at each location are independent across time.
  lbtests_gls_pval = matrix(0, nrow = length(lon), ncol = length(lat))
  #ecdfs = matrix(0, nrow = 99, ncol = n.season)
  ecdf.prop = seq(0.01, 0.99, len = 99)
  for(i in 1:length(lon))
  {
    for(j in 1:length(lat))
    {
      lbtests_gls_pval[i,j] = Box.test(R.gls.tilde[i,j,],type = "Ljung-Box",lag = 1)$p.value
      #ecdfs[,i] = sapply(ecdf.prop, function(y) mean(lbtests_pval[,i] <= y, na.rm = TRUE))
    }
  }

  ecdfs = sapply(ecdf.prop, function(y) mean(lbtests_gls_pval <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = 'Empirical Cumluative Distribution',
       xlab = "x", ylab = "ECDF",type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  
  lbtests_gwr_pval = matrix(0, nrow = length(lon), ncol = length(lat))
  ecdf.prop = seq(0.01, 0.99, len = 99)
  for(i in 1:length(lon))
  {
    for(j in 1:length(lat))
    {
      lbtests_gwr_pval[i,j] = Box.test(R.gwr.tilde[i,j,],type = "Ljung-Box",lag = 10)$p.value
    }
  }
  ecdfs = sapply(ecdf.prop, function(y) mean(lbtests_gwr_pval <= y, na.rm = TRUE))
  plot(ecdf.prop,ecdfs,main = 'Empirical Cumluative Distribution',
       xlab = "x", ylab = "ECDF",type="l")
  abline(0,1,col="blue")
  lines(ecdf.prop, ecdf.prop + 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  lines(ecdf.prop, ecdf.prop - 2*sqrt(ecdf.prop*(1-ecdf.prop)/n.loc),col="red")
  
  # pval_sort<-sort(lbtests_gwr_pval)
  # plot(pval_sort, (1:n.loc)/n.loc, type = 's', ylim = c(0, 1), 
  #      xlab = 'sample', ylab = '', main = 'Empirical Cumluative Distribution')
  # U2=(1:n.loc/n.loc)+1.96*sqrt( (1:n.loc/n.loc)*(1-1:n.loc/n.loc)/n.loc )
  # L2=(1:n.loc/n.loc)-1.96*sqrt( (1:n.loc/n.loc)*(1-1:n.loc/n.loc)/n.loc )
  # lines(pval_sort, L2, col="red")
  # lines(pval_sort, U2, col="red")
  
  #Compute threshold for CoPE sets.
  #Ignore values not on the land.
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
  level= 0.1
  norm_diff_gls = (beta.gls.wide - level) / se.gls.wide
  norm_diff_gwr = (beta.gwr.wide - level) / se.gwr.wide
  
  image.map(lon,lat,beta.gls.wide,mask=mask,col = tim.colors(64),
            horizontal=FALSE,ylab='',xlab='',
            main= paste0("GLS ",season.var[s]," level=",level))
  drawContour(lon,lat,beta.gls.wide,c=level,col="purple",lty=1)
  drawContour(lon,lat,norm_diff_gls,c=a_CoPE_gls,col="darkred")
  drawContour(lon,lat,norm_diff_gls,c=-a_CoPE_gls,col="darkgreen")
  
  image.map(lon,lat,beta.gwr.wide,mask=mask,col = tim.colors(64),
            horizontal=FALSE,ylab='',xlab='',
            main= paste0("GWR ",season.var[s]," level=",level))
  drawContour(lon,lat,beta.gwr.wide,c=level,col="purple",lty=1)
  drawContour(lon,lat,norm_diff_gwr,c=a_CoPE_gwr,col="darkred")
  drawContour(lon,lat,norm_diff_gwr,c=-a_CoPE_gwr,col="darkgreen")
}





# Compute contour of estimate. (French 2017)
x.grid = seq(0,1, length = length(lon))
y.grid = seq(0,1, length = length(lat))
cont.gls <- contourLines(x.grid,y.grid,beta.gls.wide,nlevels=1,level=level)
cont.gwr <- contourLines(x.grid,y.grid,beta.gwr.wide,nlevels=1,level=level)
#Compute a using Multiplier Bootstrap.
a_MB_gls = quantile(MBContour(x=x.grid, y=y.grid, R=R.gls.tilde, N=N, cont=cont.gls),
                probs=1 - alpha_CoPE, type=8)
a_MB_gwr = quantile(MBContour(x=x.grid, y=y.grid, R=R.gwr.tilde, N=N, cont=cont.gwr),
                probs=1 - alpha_CoPE, type=8)
#rej.region = get.rej2(a_MB,beta.wide, R)

#Compute quantile (bootstrap).
P_MC_gls = MC_gauss(R.gls.tilde,N)
P_MC_gwr = MC_gauss(R.gwr.tilde,N)
a_CoPE_gls = a_CoPE_gwr =0
while(P_MC_gls(a_CoPE_gls)>alpha_CoPE) a_CoPE_gls = a_CoPE_gls+0.01
while(P_MC_gwr(a_CoPE_gwr)>alpha_CoPE) a_CoPE_gwr = a_CoPE_gwr+0.01

#Compute quantile (Taylor's method).
Tay_fun = TaylorContour(x = x.grid,y = y.grid,cont=cont.gls,R=R.gls.tilde) #The two is for the absolute value.
a_Tay_true = 0
while(2*Tay_fun(a_Tay_true)>alpha_CoPE) a_Tay_true = a_Tay_true + 0.01
## Plot from pacakge

fields::image.plot(x.grid,y.grid,beta.gls.wide)
DrawContour(x.grid,y.grid,beta.gls.wide,level=level,col="purple", lty = 2)
DrawContour(x.grid,y.grid,norm_diff_gls,level=a_MB_gls,col="darkred")
DrawContour(x.grid,y.grid,z.wide,level=-a_MB,col="darkgreen")

n = 30
Data = ToyNoise1(n = n)
Data$z = Data$z + rep(ToySignal()$z, n)
CopeSet = ComputeCope(Data,level=4/3, mu=ToySignal()$z)
PlotCope(CopeSet)




#Apply MC simulations.
P_MC = MC_gauss(R_tilde,N)

#Compute quantile.
a_CoPE=0
while(P_MC(a_CoPE)>alpha_CoPE) a_CoPE = a_CoPE+0.05

#Ignore values not on land.
beta_hat[,,1][mask==0] = -10
#Plot the CoPE and FARE sets.
####
#Setting level.
c=0.2
#Normalized function.
norm_diff = sqrt(n)/2 * (beta_hat[,,1] - c) / sigma_hat 

setEPS()
#Plot CoPE sets.
postscript("CoPE_ks_exp.eps")
par(mar=c(3,3,2,2)+0.1)
image.map(lon,lat,beta_hat[,,1],
          mask=mask,col = tim.colors(64),horizontal=FALSE,ylab='',xlab='')
drawContour(lon,lat,beta_hat[,,1],c=c,col="purple",lty=1)
drawContour(lon,lat,norm_diff,c=a_CoPE,col="darkred")
drawContour(lon,lat,norm_diff,c=-a_CoPE,col="darkgreen")
dev.off()