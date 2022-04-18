## NA-CORDEX Util.R

## Durbin-Watson test statistic using residuals
durbinWatsonTest.r <- function(resid, max.lag=1){
  n <-  length(resid)
  dw <- rep(0, max.lag)
  den <- sum(resid^2)
  for (lag in 1:max.lag){
    dw[lag] <- (sum((resid[(lag+1):n] - resid[1:(n-lag)])^2))/den
  }
  dw
}

## Durbin-Watson test statistic using object from glm
durbinWatsonTest.glm <- function(model, max.lag=1, simulate=TRUE, reps=1000, 
                                 method="resample", 
                                 alternative="two.sided"){
  resid <- residuals(model)
  if (any(is.na(resid))) stop ('residuals include missing values')
  n <- length(resid)
  r <- dw <-rep(0, max.lag)
  den <- sum(resid^2)
  for (lag in 1:max.lag){
    dw[lag] <- (sum((resid[(lag+1):n] - resid[1:(n-lag)])^2))/den
    r[lag] <- (sum(resid[(lag+1):n]*resid[1:(n-lag)]))/den
  }
  if (!simulate){
    result <- list(r=r, dw=dw)
    result
  }
  else {
    X <- model.matrix(model)
    mu <- fitted.values(model)
    Y <-  matrix(sample(resid, n*reps, replace=TRUE), n, reps) + matrix(mu, n, reps)
    E <- apply(Y, 2, function(y) residuals(glm(y[y>0] ~ X[y>0,], family=Gamma(link = "log"))))
    if(is.list(E))
    {
      DW = unlist(lapply(E, function(x) durbinWatsonTest.r(x,max.lag=max.lag)))
    }else
      DW <- apply(E, 2, durbinWatsonTest.r, max.lag=max.lag)
    if (max.lag == 1) DW <- rbind(DW)
    p <- rep(0, max.lag)
    if (alternative == 'two.sided'){
      for (lag in 1:max.lag) {
        p[lag] <- (sum(dw[lag] < DW[lag,]))/reps
        p[lag] <- 2*(min(p[lag], 1 - p[lag]))
      }
    }
    else if (alternative == 'positive'){
      for (lag in 1:max.lag) {
        p[lag] <- (sum(dw[lag] > DW[lag,]))/reps
      }
    }
    else {
      for (lag in 1:max.lag) {
        p[lag] <- (sum(dw[lag] < DW[lag,]))/reps
      }
    }
    result <- list(r=r, dw=dw, p=p, alternative=alternative)
    result
  }
}

## Plot lines in a data list
plot_list = function(x, main = "", xlab = ""){
  maxy = sapply(x, function(j) j$y)
  plot(x[[1]], ylim = range(maxy), xlab = xlab, main = main, xlim = c(-3,3),
       col = "darkgray", lwd = 0.5)
  for (j in seq_along(x)) {
    lines(x[[j]], col = "darkgray", lwd = 0.5)
  }
}

## Multiplier Bootstrap in CoPE approach (R package CoPE)
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

## Q-Q plot with confidence band
qqnorm_new = function(y, pch = 16, xlab = "Standard Normal Quantiles",
                      add = FALSE, cex = 0.6, conf = 0.95,
                      xlim = NULL, ylim=NULL,
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
    if(is.null(xlim)) 
      xlim = range(lower-1, z-1, z+1, upper+1, na.rm = TRUE)
    if(is.null(ylim))
      ylim = xlim
    plot(z, y, type="n", xlab = xlab, ylab = ylab, pch = pch, cex = cex,
         xlim = xlim, ylim = ylim, ...)
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

## Simultaneous drawing CoPE sets at varying levels.
Simul_Cope = function(map, par, par.se, resids, level.set, zlim,
                      col.map, sub.loc, select.ind, type = "temp", 
                      alpha_CoPE = 0.1, alpha_FARE = 0.05,  N=1000)
{
  lon = unique(map$longitude)
  lon = lon[order(lon)]
  lat = unique(map$latitude)
  x.grid = seq(0,1, length = length(lon))
  y.grid = seq(0,1, length = length(lat))
  n = dim(resids)[1]

  n.levels = length(level.set)
  n.season = 4
  season.var = c("winter","spring","summer","fall" )
  for(l in 1:n.levels)
  {
    par(mfrow = c(2,2), oma = c(0, 0, 0, 3), mar = c(2, 1, 1, 2))
    for(s in 1:n.season)
    {
      
      beta.long = data.frame(lon = map$longitude, 
                             lat = map$latitude,
                             b.gwr = par[,s]
      )
      se.long = data.frame(lon = map$longitude, 
                           lat = map$latitude,
                           b.gwr = par.se[,s]
      )
      
      beta.gwr.wide = acast(beta.long, lon~lat, value.var='b.gwr')
      se.gwr.wide = acast(se.long, lon~lat, value.var='b.gwr')
      mask = (!is.na(beta.gwr.wide))*1
      
      #Compute the residuals.
      R.gwr = array(0,c(length(lon),length(lat),n))
      for(i in 1:n)
      {
        long.t = data.frame(lon = map$longitude, 
                            lat = map$latitude,
                            R.gwr = resids[i,s,])
        R.gwr[,,i] = acast(long.t, lon~lat, value.var='R.gwr')
      }
      
      sigma.gwr.hat = apply(R.gwr,1:2,function(x) sd(x, na.rm = TRUE))
      R.gwr.tilde = R.gwr / rep(sigma.gwr.hat,n)
      
      sigma.gwr.long = apply(resids[,s,],2,function(x) sd(x, na.rm = TRUE))
      R.gwr.long = resids[,s,]/sigma.gwr.long
      R.gwr.sub = R.gwr.long[,select.ind]
      
      A = (mask==1)
      A[is.na(A)] = FALSE
      for(i in 1:n) R.gwr.tilde[,,i] = R.gwr.tilde[,,i] * A
      
      #Compute quantile (bootstrap).
      P_MC_gwr = MC_gauss(R.gwr.tilde,N)
      a_CoPE_gwr = 0
      while(P_MC_gwr(a_CoPE_gwr)>alpha_CoPE/n.levels) a_CoPE_gwr = a_CoPE_gwr+0.01
      
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
      
      gwr_interp = mba.surf(cbind(map[,c("longitude","latitude")],
                                  norm_diff_gwr), 200, 200)$xyz.est
      #zlim = c(0,0.6)
      slope_interp <- mba.surf(cbind(map[,c("longitude","latitude")],
                                    par[,s]), 200, 200)$xyz.est
      mask_large = (!is.na(slope_interp$z))*1
      
      image.map(slope_interp$x,slope_interp$y,
                slope_interp$z,mask=mask_large,col = col.map,
                #horizontal=FALSE,
                legend = FALSE,
                ylab='',xlab='',
                zlim = zlim,
                main= paste0(season.var[s]))
      maps::map("county",add = TRUE)
      points(sub.loc$x, sub.loc$y,pch = 16, col = "black",cex=1.1)
      drawContour(slope_interp$x,slope_interp$y, bound = FALSE,
                  slope_interp$z,c=level,col="purple",lty=1)
      if(type=="temp")
      {
        drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=a_CoPE_gwr,col="darkred")
        drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=-a_CoPE_gwr,col="darkgreen")
      }
      if(type=="prep")
      {
        drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=a_CoPE_gwr,col="darkgreen")
        drawContour(gwr_interp$x,gwr_interp$y,gwr_interp$z,c=-a_CoPE_gwr,col="darkred")
      }
      
    }
    par(mfrow=c(1, 1), mar=c(1, 1, 0.7, 2), new=FALSE)
    image.plot(zlim= zlim, legend.only=TRUE, col=col.map, legend.mar = 0.1,
               legend.width=0.8, legend.shrink=1,legend.cex = 0.6)
  }
}

# The following functions were modified from Sommerfield et al. (2018)
## This function takes realizations Y of a random field and returns a distribution functions
## for the supremum of the limiting Gaussian field.
MC_gauss = function(Y,N){
  n = dim(Y)[3]
  Y = matrix(Y,ncol=n)
  g = matrix(rnorm(n*N),n,N)
  maxima = apply(abs(Y %*% g), 2, function(x) max(x, na.rm = TRUE)) / sqrt(n)
  
  function(t) sum(maxima>=t) / length(maxima)
}

## Plotting values on the map.
image.map = function(lon, lat, img,legend = TRUE, mask=NULL, xlab='longitude', ylab='latitude', ...) {
  if(!is.null(mask)) {
    img = img*mask
    xlim = lon[range(which(rowSums(mask, na.rm=TRUE)>0))]
    ylim = lat[range(which(colSums(mask, na.rm=TRUE)>0))]
  } else {
    xlim = range(lon)
    ylim = range(lat)
  }
  if(legend) image.plot(lon, lat, img, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  else image(lon, lat, img, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  #map('world', add=TRUE)
}

## Drawing the contour.
drawContour = function(x,y,z,c,col,lty=1, bound = TRUE, alternative="greater"){
  if(alternative == "greater")
  {
    u.loc = which(z == max(z,na.rm = TRUE), arr.ind = TRUE)
    l.loc = which(z == min(z,na.rm = TRUE), arr.ind = TRUE)
    if(bound==FALSE & all(z>c,na.rm = TRUE))
    {
      # text(gls_interp$x[u.loc[1]]-0.5, gls_interp$y[u.loc[2]]-0.5,
      #      "ES",xpd=NA, col=col, cex=1.5)
      text(x[which.min(x)]+(max(x)-min(x))/8, 
           y[which.min(y)]+2*(max(y)-min(y))/6,
           "ES-H",xpd=NA, col=col, cex=1.3) 
    }else if(bound==FALSE & all(z<c,na.rm = TRUE))
    {
      # text(gls_interp$x[u.loc[1]]-0.5, gls_interp$y[u.loc[2]]-0.5,
      #      "ES",xpd=NA, col=col, cex=1.5)
      text(x[which.min(x)]+(max(x)-min(x))/8, 
           y[which.min(y)]+2*(max(y)-min(y))/6,
           "ES-L",xpd=NA, col=col, cex=1.3) 
    }else if(c>=0 & all(z>c,na.rm = TRUE))
    {
      #print("Full Set")
      # text(gls_interp$x[u.loc[1]]-0.5, gls_interp$y[u.loc[2]]+0.5,
      #      "Full Upper Set", xpd=NA, col=col)
      text(x[which.min(x)]+(max(x)-min(x))/8, 
           y[which.min(y)]+3*(max(y)-min(y))/6,
           "FUS",xpd=NA, col=col, cex=1.3)
    }else if(c<=0 & all(z<c,na.rm = TRUE))
    {
      text(x[which.min(x)]+(max(x)-min(x))/8, 
           y[which.min(y)]+(max(y)-min(y))/6,
           #"Full Lower Set", xpd=NA, col=col
           "FLS",xpd=NA, col=col, cex=1.3)
    }else if(c>=0 & all(z<c,na.rm = TRUE))
    {
      # text(gls_interp$x[u.loc[1]]-0.5, gls_interp$y[u.loc[2]]-0.5,
      #      #"Empty Upper Set", xpd=NA, col=col
      #      "EUS",xpd=NA, col=col, cex=1.5)
      text(x[which.min(x)]+(max(x)-min(x))/8, 
           y[which.min(y)]+3*(max(y)-min(y))/6,
           "EUS",xpd=NA, col=col, cex=1.3)
    }else if(c<=0 & all(z>c,na.rm = TRUE))
    {
      #print("Empty Set")
      text(x[which.min(x)]+(max(x)-min(x))/8, 
           y[which.min(y)]+(max(y)-min(y))/6,
           "ELS",xpd=NA, col=col, cex=1.3)
    }else
    {
      C<-contourLines(x,y,z,levels=c,nlevels=1)
      
      for(i in 1:length(C))
        lines(C[[i]]$x,C[[i]]$y,pch=20,col=col,lwd=3,lty=lty)
    }
  }else if(alternative == "less")
  {
    if(c>=0 & all(z>c,na.rm = TRUE))
    {
      #print("Full Set")
      text(max(x)-1, max(y)-1, "Full Lower Set", xpd=NA, col=col)
    }else if(c<=0 & all(z<c,na.rm = TRUE))
    {
      text(min(x)+1, min(y)+1, "Full Upper Set", xpd=NA, col=col)
    }else if(c>=0 & all(z<c,na.rm = TRUE))
    {
      text(max(x)-1, max(y)-1, "Empty Lower Set", xpd=NA, col=col)
    }else if(c<=0 & all(z>c,na.rm = TRUE))
    {
      #print("Empty Set")
      text(min(x)+1, min(y)+1, "Empty Upper Set", xpd=NA, col=col)
    }else
    {
      C<-contourLines(x,y,z,levels=c,nlevels=1)
      
      for(i in 1:length(C))
        lines(C[[i]]$x,C[[i]]$y,pch=20,col=col,lwd=3,lty=lty)
    }
  }
}