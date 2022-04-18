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

mq_print <- function(x,lag=24,adj=0){
  # Compute multivariate Ljung-Box test statistics
  #
  # adj: adjustment for the degrees of freedomm in the chi-square distribution.
  # adj is the number of coefficient parameters used in the fitted model, if any.
  #
  if(!is.matrix(x))x=as.matrix(x)
  nr=nrow(x)
  nc=ncol(x)
  g0=var(x)
  ginv=solve(g0)
  qm=0.0
  QM=NULL
  df = 0
  for (i in 1:lag){
    x1=x[(i+1):nr,]
    x2=x[1:(nr-i),]
    g = cov(x1,x2)
    g = g*(nr-i-1)/(nr-1)
    h=t(g)%*%ginv%*%g%*%ginv
    qm=qm+nr*nr*sum(diag(h))/(nr-i)
    df=df+nc*nc
    dff= df-adj
    mindeg=nc^2-1
    pv = 1
    if(dff > mindeg)pv=1-pchisq(qm,dff)
    QM=rbind(QM,c(i,qm,dff,pv))
  }
  pvs=QM[,4]
  dimnames(QM) = list(names(pvs),c("  m  ","    Q(m) ","   df  "," p-value"))
  return(QM)
  # cat("Ljung-Box Statistics: ","\n")
  # printCoefmat(QM,digits = 3)
  # #
  # par(mfcol=c(1,1))
  # plot(pvs,ylim=c(0,1),xlab="m",ylab="prob",main="p-values of Ljung-Box statistics")
  # abline(h=c(0))
  # lines(rep(0.05,lag),lty=2,col='blue')
}

varxfixed_pred <- function(m1,newxt=NULL,hstep=1,orig=0){
  #This program predicts the VARX model.
  ## z(t) = c0 + sum_{i=1}^p phi_i * z(t-i) + \sum_{j=0}^m xt(t-j) + a(t).
  ##
  zt=as.matrix(m1$data); xt=as.matrix(m1$xt); p=m1$aror; m=m1$m
  Ph0=as.matrix(m1$Ph0); Phi=as.matrix(m1$Phi); Sig=as.matrix(m1$Sigma); beta=as.matrix(m1$beta)
  include.m=m1$include.mean
  nT=dim(zt)[1]; k=dim(zt)[2]; dx=dim(xt)[2]
  se=NULL
  if(length(Ph0) < 1)Ph0=matrix(rep(0,k),k,1)
  if(hstep < 1)hstep=1
  if(orig < 1) orig=nT
  #
  if(length(newxt) > 0){
    if(!is.matrix(newxt))newxt=as.matrix(newxt)
    ### calculate predictions
    h1=dim(newxt)[1]
    hstep=min(h1,hstep)
    nzt=as.matrix(zt[1:orig,])
    if(dx > 1){
      xt=rbind(xt[1:orig,,drop=FALSE],newxt)
    }else{
      xt=as.matrix(c(c(xt[1:orig,]),c(newxt)))
    }
    for (i in 1:hstep){
      tmp=Ph0
      ti=orig+i
      ## VAR part
      for (i in 1:p){
        idx=(i-1)*k
        tmp=tmp+Phi[,(idx+1):(idx+k)]%*%matrix(nzt[ti-i,],k,1)
      }
      if(m > -1){
        for (j in 0:m){
          jdx=j*dx
          tmp=tmp+beta[,(jdx+1):(jdx+dx)]%*%matrix(xt[ti-j,],dx,1)
        }
      }
      nzt=rbind(nzt,c(tmp))
    }
    ### compute standard errors of predictions
    mm=VARpsi(Phi,lag=hstep)
    Si=Sig
    se=matrix(sqrt(diag(Si)),1,k)
    if(hstep > 1){
      for (i in 2:hstep){
        idx=(i-1)*k
        wk=as.matrix(mm$psi[,(idx+1):(idx+k)])
        Si=Si+wk%*%Sig%*%t(wk)
        se1=sqrt(diag(Si))
        se=rbind(se,se1)
      }
    }
  }else{
    cat("Need new data for input variables!","\n")
  }
  #
  return( list(pred=nzt[(orig+1):(orig+hstep),],se=se[1:hstep,],orig=orig,h=hstep))
}
