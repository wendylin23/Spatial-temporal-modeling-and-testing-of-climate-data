library(matlib)

varx_fixed = function(zt,p,xt,fixed=NULL,m=0,include.mean=T)
{
  zt=as.matrix(zt)
  if(length(xt) < 1){
    m = -1; kx=0}else{
      xt=as.matrix(xt); kx=dim(xt)[2]
    }
  if(p < 0)p=0
  ist=max(p,m)+1
  nT=dim(zt)[1]
  k=dim(zt)[2]
  yt=zt[ist:nT,]
  xmtx=NULL
  if(include.mean)xmtx=rep(1,(nT-ist+1))
  #
  if(p > 0){
    for (i in 1:p){
      xmtx=cbind(xmtx,zt[(ist-i):(nT-i),])
    }
  }
  #
  if( m > -1){
    for (j in 0:m){
      xmtx=cbind(xmtx,xt[(ist-j):(nT-j),])
    }
  }
  p1=dim(xmtx)[2]
  nobe=dim(xmtx)[1]
  if(length(fixed) < 1){
    ## no constriants
    xpx=t(xmtx)%*%xmtx
    xpy=t(xmtx)%*%yt
    xpxi=solve(xpx)
    beta=xpxi%*%xpy
    resi=as.matrix(yt-xmtx%*%beta)
    sig=crossprod(resi,resi)/nobe
    co=kronecker(sig,xpxi)
    se=sqrt(diag(co))
    se.beta=matrix(se,nrow(beta),k)
    npar=nrow(beta)*k
    d1=log(det(sig))
    aic=d1+2*npar/nobe
    bic=d1+(npar*log(nobe))/nobe
  }else{
    beta=matrix(0,p1,k)
    se.beta=matrix(1,p1,k)
    resi=yt
    npar=0
    coef.cov = list()
    for (i in 1:k){
      idx=c(1:p1)[fixed[,i] > 0]
      npar=npar+length(idx)
      if(length(idx) > 0){
        xm=as.matrix(xmtx[,idx])
        y1=matrix(yt[,i],nobe,1)
        xpx=t(xm)%*%xm
        xpy=t(xm)%*%y1
        #xpxi=solve(xpx)
        xpxi=Ginv(xpx)
        beta1=xpxi%*%xpy
        res = y1 - xm%*%beta1
        sig1=sum(res^2)/nobe
        se=sqrt(diag(xpxi)*sig1)
        beta[idx,i]=beta1
        se.beta[idx,i]=se
        resi[,i]=res
        coef.cov[[i]]= (xpxi * sig1)
      }
    }
    sig=crossprod(resi,resi)/nobe
    d1=log(det(sig))
    aic=d1+2*npar/nobe
    bic=d1+log(nobe)*npar/nobe
    # end of for (i in 1:k)
  }
  
  Ph0=NULL
  icnt=0
  if(include.mean){
    Ph0=beta[1,]; icnt=icnt+1
  }
  Phi=NULL
  sePhi=NULL
  if(p > 0){
    Phi=t(beta[(icnt+1):(icnt+k*p),])
    sePhi=t(se.beta[(icnt+1):(icnt+k*p),])
    icnt=icnt+k*p
    ## end of if(p > 0)
  }
  if(m > -1){
    Beta=t(beta[(icnt+1):(icnt+(m+1)*kx),])
    seBeta=t(se.beta[(icnt+1):(icnt+(m+1)*kx),])
    if(kx == 1){
      Beta=t(Beta)
      seBeta=t(seBeta)
    }
    ## end of if(m > -1)
  }
  return(list(data=zt,xt=xt,aror=p,m=m,
              Ph0=Ph0,Phi=Phi,sePhi = sePhi,
              beta=Beta,se.beta = seBeta,
              coef=beta,se.coef=se.beta,coef.cov = coef.cov,
              residuals=resi,Cor = cov2cor(sig),Sigma=sig,
              aic = aic,bic = bic,
              include.mean=include.mean))
}

