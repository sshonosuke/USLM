#=======================================================================
#  R code for 
#   "Small Area Estimation via Unmatched Sampling and Linking Models"
#
#       By Shonosuke Sugasawa, Tatsuya Kubokawa and J. N. K. Rao
#=======================================================================

# ABBREVIATIONS
# USLM: unmatched sampling and linking model
# GH: Gauss-Hermite 
# MSE: mean squared error
# MH: Metropolis-Hastings


#-----------------------------------------------------------------------
#  Functions for fitting parametric USLM 
#-----------------------------------------------------------------------

# Function: USLM
#
# INPUT:
# y: A vector of responce variables.
# X: A matrix of covariates whose column vectors have the same dimension as y.
# Di: A vector of sampling variances.
# link: A link function chosen from "exp" (exponetial link) or "logit" (logit link). 
#       The default is "exp".
# n: The number of GH quadrature. The default value is 100.
# dd: A small positive value. The EM algorithm is terminated if the relative 
#     difference is smaller than this value. Default value is 0.0001
# Init: A vector of initial parameter values. The default is NA in which the 
#     initial vlaues are set to A=1 and 0 for regression coeffieicnts.
#
# OUTPUT: 
# Est: A vector of estimates of model parameters.
# EB: A vector of empirical Bayes estimates.
# AEB: A vector of approximate empirical Bayes estimates.


USLM=function(y,X,Di,link="exp",n=100,dd=0.0001,Init=NA){
  library(statmod)
  m=length(y); p=dim(X)[2]
  if(link=="exp"){ Ln=function(x){ exp(x) } }
  if(link=="logit"){ Ln=function(x){ exp(x)/(1+exp(x)) } }
  if(is.na(Init[1])){ Init=c(rep(0,p),1) }
  
  Pred=function(para){
    rr=gauss.quad(n,kind="hermite")
    ww=rr$weights/sqrt(pi); u=rr$nodes
    beta=para[1:p]; A=para[p+1]; mu=as.vector(X%*%beta)
    z=t(matrix(rep(u,m),n,m)); U=sqrt(2)*sqrt(A)*z+mu
    th=Ln(U); Den=t(exp(-(2*Di)^(-1)*(y-th)^2))
    pv=apply(t(th)*Den*ww,2,sum)/apply(Den*ww,2,sum)
    return(pv)
  }
  
  Approx.Pred=function(para){
    beta=para[1:p]; A=para[p+1]; mu=as.vector(X%*%beta)
    v.ast=c()
    for(i in 1:m){
      pos=function(v){ v^2/A+(y[i]-Ln(mu[i]+v))^2/Di[i] }
      v.ast[i]=optimize(f=pos,interval=c(-100,100))$minimum
    }
    return(Ln(mu+v.ast))
  }
  
  E.step=function(para){
    rr=gauss.quad(n,kind="hermite")
    ww=rr$weights/sqrt(pi); u=rr$nodes
    beta=para[1:p]; A=para[p+1]; mu=as.vector(X%*%beta)
    z=t(matrix(rep(u,m),n,m)); U=sqrt(2)*sqrt(A)*z+mu
    th=Ln(U); Den=t(exp(-(2*Di)^(-1)*(y-th)^2))
    I1=apply(t(U)*Den*ww,2,sum)/apply(Den*ww,2,sum)
    I2=apply(t(U^2)*Den*ww,2,sum)/apply(Den*ww,2,sum)
    return(cbind(I1,I2))
  }
  
  beta=Init[1:p]; A=Init[p+1]
  d=1
  while(d>dd){
    Int=E.step(c(beta,A))
    int1=Int[,1]; int2=Int[,2]
    new.beta=solve(t(X)%*%X)%*%t(X)%*%int1
    new.A=mean(int2-2*X%*%new.beta*int1+(X%*%new.beta)^2)
    d=sum(abs(c(new.beta-beta,new.A-A)))/(sum(abs(c(beta,A)))+0.0001)
    beta=new.beta; A=max(new.A,0)
  }
  
  Est=c(beta,A); names(Est)=c(paste0("beta",0:(p-1)),"A")
  EB=Pred(Est); AEB=Approx.Pred(Est)
  Result=list(Est,EB,AEB); names(Result)=c("Est","EB","AEB")
  return(Result)
}






#-----------------------------------------------------------------------
#  Functions for estimating area-specific MSE in parametric USLM 
#-----------------------------------------------------------------------

# Function: mseUSLM
#
# INPUT:
# y: A vector of responce variables.
# X: A matrix of vovariates whose column vectors have the same dimension as y.
# Di: A vector of sampling variances.
# link: A link function chosen from "exp" (exponetial link) or "logit" (logit link). 
#     The default is "exp".
# n: The number of GH quadrature. The default value is 100.
# dd: A small positive value. The EM algorithm is terminated if the relative 
#     difference is smaller than this value. Default value is 0.0001
# Init: A vector of initial parameter values. The default is NA in which the 
#     initial vlaues are set to A=1 and 0 for regression coeffieicnts.
# mc: The number of Monte Carlo samples used for computing the leading term of MSE.
#     The default value is 1000
# B: The number of bootstrap samples for MSE
#
# OUTPUT: A vector of MSE estimates.

mseUSLM=function(y,X,Di,link="exp",n=100,dd=0.0001,Init=NA,mc=1000,B=100){
  library(statmod)
  rr=gauss.quad(n,kind="hermite")
  ww=rr$weights/sqrt(pi); u=rr$nodes
  
  m=length(y); p=dim(X)[2]
  if(link=="exp"){ Ln=function(x){ exp(x) } }
  if(link=="logit"){ Ln=function(x){ exp(x)/(1+exp(x)) } }
  if(is.na(Init[1])){ Init=c(rep(0,p),1) }
  
  pred1=function(y,para){
    beta=para[1:p]; A=para[p+1]; mu=as.vector(X%*%beta)
    z=t(matrix(rep(u,m),n,m)); th=Ln(sqrt(2)*sqrt(A)*z+mu)
    Den=t(exp(-(2*Di)^(-1)*(y-th)^2))
    pv=apply(t(th)*Den*ww,2,sum)/apply(Den*ww,2,sum)
    return(pv)
  }
  
  pred2=function(yy,num,para){
    len=length(yy)
    beta=para[1:p]; A=para[p+1]
    mu=(X%*%beta)[num]; DD=Di[num]
    z=t(matrix(rep(u,len),n,len)); th=Ln(sqrt(2)*sqrt(A)*z+mu)
    Den=t(exp(-(2*DD)^(-1)*(yy-th)^2))
    pv=apply(t(th)*Den*ww,2,sum)/apply(Den*ww,2,sum)
    return(pv)
  }
  
  T1=function(para){
    beta=para[1:p]; A=para[p+1]; mu=X%*%beta
    val=c()
    for(i in 1:m){
      th=Ln(sqrt(A)*rnorm(mc)+mu[i])
      yy=th+rnorm(mc,0,sqrt(Di[i]))
      BP=pred2(yy,i,para)
      val[i]=mean((BP-th)^2)
    }
    return(val)
  }
  
  est=USLM(y,X,Di,link=link,n=n,dd=dd,Init=Init)[[1]]
  hbeta=est[1:p]; hA=est[p+1]
  hmu=as.vector(X%*%hbeta)
  T1.boot=matrix(NA,B,m); T2.boot=matrix(NA,B,m)
  for(b in 1:B){
    yb=Ln(hmu+rnorm(m,0,sqrt(hA)))+rnorm(m,0,sqrt(Di))
    boot.est=USLM(yb,X,Di,link=link,n=n,dd=dd,Init=Init)[[1]]
    T2.boot[b,]=(pred1(yb,boot.est)-pred1(yb,est))^2
    T1.boot[b,]=T1(boot.est)
  }
  
  MSE=2*T1(est)-apply(na.omit(T1.boot),2,mean)+apply(na.omit(T2.boot),2,mean)
  MSE=ifelse(MSE>0,1,0)*MSE
  return(MSE)
}




#-----------------------------------------------------------------------
#  Functions for fitting USLM with P-spline 
#-----------------------------------------------------------------------

# Function: spUSLM
#
# INPUT:
# y: A vector of responce variables.
# X: A matrix of vovariates whose column vectors have the same dimension as y.
# Di: A vector of sampling variances.
# link: A link function chosen from "exp" (exponetial link) or "logit" (logit link). 
#       The default is "exp".
# dd: A small positive value. The EM algorithm is terminated if the relative 
#     difference is smaller than this value. Default value is 0.0001
# Init: A vector of initial parameter values. The default is NA in which the 
#     initial vlaues are set to A=1 and 0 for regression coeffieicnts.
# mc: The initial number of Monte Carlo samples used in E-step. The default value is 3000.
# burn: The number of burn-in period in E-step. The default value is 2000.
# band: The value of a stap-size in MH update in E-step. The default value is 0.3.
#
# OUTPUT: 
# Est: A vector of estimates of model parameters.
# EB: A vector of empirical Bayes estimates.
# AEB: A vector of approximate empirical Bayes estimates.
# Gamma: A vector of estimates of gamma (coefficients of spline terms).


spUSLM=function(y,X,Z,Di,link="exp",dd=0.005,Init=NA,mc=3000,burn=2000,band=0.3){
  library(MASS)
  library(MCMCpack)
  m=length(y); p=dim(X)[2]; K=dim(Z)[2]
  if(link=="exp"){ Ln=function(x){ exp(x) } }
  if(link=="logit"){ Ln=function(x){ exp(x)/(1+exp(x)) } }
  if(is.na(Init[1])){ Init=c(rep(0,p),1,1) }
  
  approx=function(para){
    beta=para[1:p]; A=para[p+1]; tau=para[p+2]
    mu=as.vector(X%*%beta)
    v.ast=rep(0,m); gam.ast=rep(0,K)
    d=1; count=1
    while(d>dd){
      old=v.ast
      mu2=as.vector(Z%*%gam.ast)
      for(i in 1:m){
        pos.v=function(a){ a^2/A+(y[i]-Ln(mu[i]+mu2[i]+a))^2/Di[i] }
        v.ast[i]=optim(par=0,fn=pos.v,method="BFGS")$par
      }
      pos.gam=function(a){ sum(a^2)/tau^2+sum((y-Ln(mu+Z%*%a+v.ast))^2/Di) }
      gam.ast=optim(par=gam.ast,fn=pos.gam,method="BFGS")$par
      d=sum(abs(v.ast-old))/sum(abs(old))
    }
    return(list(v.ast,gam.ast))
  }
  
  logdens.u=function(u,beta,A,gam){
    val=-0.5*(y-Ln(u))^2/Di-0.5*(u-X%*%beta-Z%*%gam)^2/A
    return(as.vector(val))
  }
  
  tr=function(mat){ sum(diag(mat)) }
  
  beta=Init[1:p]; A=Init[p+1]; tau=Init[p+2]
  d=1; count=1
  while(d>dd){
    if(count>30){ mc=mc+500 }
    gam.pos=matrix(NA,mc,K); u.pos=matrix(NA,mc,m); e.pos=matrix(NA,mc,m)
    gam.pos[1,]=rep(0,K); u.pos[1,]=rep(0,m)
    mat=solve(t(Z)%*%Z/A+diag(K)/tau^2)
    for(k in 2:mc){
      uu=u.pos[k-1,]
      gam.pos[k,]=mvrnorm(1,mat%*%t(Z)%*%(uu-X%*%beta)/A,mat)
      ggam=gam.pos[k,]
      prop=uu+rnorm(m,0,band)
      prob1=logdens.u(uu,beta,A,ggam); prob2=logdens.u(prop,beta,A,ggam)
      rr=exp(prob2-prob1)
      for(i in 1:m){ rr[i]=min(1,rr[i]) }
      ch=rbinom(m,1,rr)
      u.pos[k,]=uu+ch*(prop-uu)
      e.pos[k,]=u.pos[k,]-Z%*%gam.pos[k,]
    }
    om=(1:burn)
    uu=u.pos[-om,]; gg=gam.pos[-om,]; ee=e.pos[-om,]
    
    new.beta=solve(t(X)%*%X)%*%t(X)%*%apply(ee,2,mean)
    hmu=as.vector(X%*%new.beta)
    
    V1=c(); V2=c()
    for(k in 2:mc){ 
      V1[k]=sum((e.pos[k,]-hmu)^2)
      V2[k]=sum((gam.pos[k,])^2)
    }
    
    v1=mean(V1[-om]); v2=mean(V2[-om])
    new.A=v1/m; new.tau=sqrt(v2/K)
    d=sum(abs(c(new.beta-beta,new.A-A,new.tau-tau)))/(sum(abs(c(beta,A,tau)))+0.0001)
    beta=new.beta; A=new.A; tau=new.tau
    count=count+1
  }
  Est=c(beta,A,tau); names(Est)=c(paste0("beta",0:(p-1)),"A","tau")
  
  gam.pos=matrix(NA,mc,K); u.pos=matrix(NA,mc,m)
  gam.pos[1,]=rep(0,K); u.pos[1,]=rep(0,m)
  mat=solve(t(Z)%*%Z/A+diag(K)/tau^2)
  
  for(k in 2:mc){
    uu=u.pos[k-1,]
    ggam=mvrnorm(1,mat%*%t(Z)%*%(uu-X%*%beta)/A,mat)
    gam.pos[k,]=ggam
    prop=uu+rnorm(m,0,band)
    prob1=logdens.u(uu,beta,A,ggam); prob2=logdens.u(prop,beta,A,ggam)
    rr=exp(prob2-prob1)
    for(i in 1:m){ rr[i]=min(1,rr[i]) }
    ch=rbinom(m,1,rr)
    u.pos[k,]=uu+ch*(prop-uu)
  }
  uu=u.pos[-om,]; gg=gam.pos[-om,]
  
  EB=apply(Ln(uu),2,mean)
  gamh=apply(gg,2,mean)
  res=approx(Est)
  AEB=as.vector(Ln(X%*%beta+Z%*%res[[2]]+res[[1]]))
  Result=list(Est,EB,AEB,gamh); names(Result)=c("Est","EB","AEB","Gamma")
  return(Result)
}



