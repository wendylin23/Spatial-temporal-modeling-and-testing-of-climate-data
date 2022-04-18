library(matlib)

## weight functions

##### from spgwr pacakge
gwr.Gauss <- function (dist2, bandwidth) {
  w <- exp((-0.5)*((dist2)/(bandwidth^2)))
  w
}

gwr.gauss <- function(dist2, bandwidth) {
  w <- exp((-dist2)/(bandwidth^2))
  w
}

gwr.bisquare <- function(dist2, d) {
  d2 <- d^2
  w <- ifelse(dist2 > d2, 0, (1 - (dist2/d2))^2)
  w
}

gwr.lm <- function(dist, d) {
  #d2 <- d^2
  w <- ifelse(dist > d, 0, (1 - (dist/d)))
  w
}

gwr.exp <- function(dist,bandwidth) {
  #d2 <- d^2
  w <- exp((-dist)/(bandwidth))
  w
}

gwr.inv <- function(dist)
{
  w <- diag(1, dim(dist)[1])
  w[dist==1] <- 0.25
  w
}

varw_fixed = function(zt,p,xt,w = NULL,fixed=NULL,m=0,include.mean=T)
{
  if(is.null(w))
  {
    zt=as.matrix(zt)
  }
  # if(length(xt) < 1){
  #   m = -1; kx=0}else{
  #     xt=as.matrix(xt); kx=dim(xt)[2]
  #   }
  if(p < 0)p=0
  ist=max(p,m)+1
  nT=dim(zt)[1]
  J=dim(zt)[2]
  n=dim(zt)[3]
  nj = n*nT
  #yt=zt[ist:nT,,]
  yt = zt 
  xmtx=xt
  #if(include.mean)xmtx=rep(1,(nT-ist+1))
  #
  if(p > 0){
    print("not ready")
    # for (i in 1:p){
    #   xmtx=cbind(xmtx,zt[(ist-i):(nT-i),])
    # }
  }
  #
  # if( m > -1){
  #   for (j in 0:m){
  #     xmtx=cbind(xmtx,xt[(ist-j):(nT-j),])
  #   }
  # }
  p1=dim(xmtx)[2]
  nobe=dim(xmtx)[1]
  if(length(fixed) < 1){
    print("not ready")
    ## no constriants
    # xpx=t(xmtx)%*%xmtx
    # xpy=t(xmtx)%*%yt
    # xpxi=solve(xpx)
    # beta=xpxi%*%xpy
    # resi=as.matrix(yt-xmtx%*%beta)
    # sig=crossprod(resi,resi)/nobe
    # co=kronecker(sig,xpxi)
    # se=sqrt(diag(co))
    # se.beta=matrix(se,nrow(beta),J)
    # npar=nrow(beta)*J
    # d1=log(det(sig))
    # aic=d1+2*npar/nobe
    # bic=d1+(npar*log(nobe))/nobe
  }else{
    beta=array(0,dim=c(p1,J,n))
    se.beta = se.beta.sw=array(1,dim=c(p1,J,n))
    resi=yt
    npar=0
    sig = array(0,dim = c(J,n,n))
    aic = rep(0,J)
    rss = rep(0,J)
    for (j in 1:J){
      idx=c(1:p1)[fixed[,j] > 0]
      npar=npar+length(idx) ##？
      #si = matrix(0,nrow=n*nT,ncol=n*nT)
      si = list()
      sig.t = rep(0, n)
      se.w = array(0, dim = c(length(idx),n))
      for(i in 1:n)
      {
        wi = w[i,]
        if(length(idx) > 0){
          xm = xmtx[,idx,]
          xi = matrix(0,nrow = length(idx),ncol = length(idx))
          xyi = matrix(0, nrow = length(idx), ncol = 1)
          #xw = NULL
          xw = list()
          for(ii in 1:n)
          {
            xi = xi + wi[ii]*t(xm[,,ii])%*%xm[,,ii]
            xyi = xyi + wi[ii]*t(xm[,,ii])%*%matrix(yt[,j,ii],nobe,1)
            #xw = cbind(xw,wi[ii]*t(xm[,,ii]))
            xw[[ii]] = wi[ii]*t(xm[,,ii])
          }
          xw = do.call("cbind", xw)
          #xpx=t(xm)%*%xm
          #xpy=t(xm)%*%y1
          #xpxi=solve(xpx)
          xxi=Ginv(xi)
          beta1=xxi%*%xyi
          res = yt[,j,i] - xm[,,i]%*%beta1
          sig1=sum(res^2)#/nobe
          sig.t[i] = sig1
          xw = xxi%*%xw
          se=sqrt(diag(xw%*%t(xw)*sig1))
          se.w[,i] = diag(xw%*%t(xw))
          #se=sqrt(diag(xxi)*sig1)
          #se=sqrt(diag(Ginv(xw%*%t(xw)))*sig1)
          beta[idx,j,i]=beta1
          se.beta[idx,j,i]=se
          resi[,j,i]=res
        }else{
          print("Please specify the correct constraints")
        }
        si[[i]] = xm[,,i]%*%xw#/nobe
      }
      si = do.call("rbind",si)
      sig.w = sum(sig.t)/(nj-tr(si))
      se.beta.sw[idx,j,] = sqrt(se.w*sig.w)
      sig[j,,] = crossprod(resi[,j,],resi[,j,])/nobe
      #d1 = log(sum(diag(sig[j,,]))/n)
      d1 = sum(resi[,j,]^2)/nj
      #d1 = log(det(sig[j,,]))
      #aic[j]=2*n*d1+ (n+tr(si))/(n+2-tr(si))#2*npar/nobe
      aic[j]=2*nj*d1+ nj * (nj+tr(si))/(nj-2-tr(si))#2*npar/nobe
      #aic[j]=2*nj*d1+ (nj+tr(si))/(nj+2-tr(si))
      rss[j] = d1
    }
    # sig = array(0,dim = c(J,J,n.loc))
    # aic = bic = rep(0, n.loc)
    # for(i in 1:n.loc)
    # {
    #   sig[,,i]=crossprod(resi[,,i],resi[,,i])/nobe
    #   d1=log(det(sig[,,i]))
    #   aic[i]=d1+2*npar/nobe
    #   #bic[i]=d1+log(nobe)*npar/nobe
    # }
    # end of for (i in 1:J)
  }
  
  # Ph0=NULL
  # icnt=0
  # if(include.mean){
  #   Ph0=beta[1,]; icnt=icnt+1
  # }
  # Phi=NULL
  # sePhi=NULL
  # if(p > 0){
  #   Phi=t(beta[(icnt+1):(icnt+J*p),])
  #   sePhi=t(se.beta[(icnt+1):(icnt+J*p),])
  #   icnt=icnt+J*p
  #   ## end of if(p > 0)
  # }
  # if(m > -1){
  #   Beta=t(beta[(icnt+1):(icnt+(m+1)*kx),])
  #   seBeta=t(se.beta[(icnt+1):(icnt+(m+1)*kx),])
  #   if(kx == 1){
  #     Beta=t(Beta)
  #     seBeta=t(seBeta)
  #   }
  #   ## end of if(m > -1)
  # }
  return(list(data=zt,xt=xt,aror=p,m=m,
              #Ph0=Ph0,Phi=Phi,sePhi = sePhi,
              #beta=Beta,se.beta = seBeta,
              coef=beta,se.coef=se.beta.sw,
              residuals=resi,
              #Cor = cov2cor(sig),
              Sigma=sig,
              aic = aic,
              rss = rss,
              #bic = bic,
              include.mean=include.mean))
}