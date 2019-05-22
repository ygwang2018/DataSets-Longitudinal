

##################################AR(1) matrix
ar1=function(rho,n)
{
  e=1:n
  rho^abs(outer(e,e,"-"))
  
}

#####################################
# EXchangeable matrix
##################################
cs<- function(rho,n) {
  lags<-matrix(1,n,n)-diag(1,n,n)
  rho^lags 
}


####################################
#MA(1) matrix
####################################
ma1=function(rho,n){
  out=diag(n)/2
  out[row(out)==col(out)+1]=rho
  out+t(out)
}


##################################################
# Toeplitz matrix
####################################################
top=function(rho){
  
  a=c(1,rho)
  toeplitz(a)
  
}
####################
#inverse matrix of exchangeable matrix
#####################################
cs.inv=function(rho,n){
  
  inv=1/(1-rho)*diag(n)-rho/((1-rho)*(1+(n-1)*rho))*matrix(1,nc=n,nr=n)
  
  return(inv)
}

####################
#inverse matrix of AR(1) matrix
#####################################
ar.inv=function(rho,n){
  
  inv=1/(1-rho^2)*(ma1(-rho,n)+diag(c(0,rep(rho^2,n-2),0)))
  return(inv) 
}



######################################################

######################################geeglm model selection criteria for QIC CIC RJ Delta1 Delta2
QIC.CIC.geeglm <- function(datatype,model.geeglm, model.independence)
  
{
  
  #calculates binomial QIC of Pan (2001)
  
  #obtain trace term of QIC~CIC
  ############################inverse matrix of covariance matrix for parameter estimation by independence working model
  
  OmegaI <- solve(model.independence$geese$vbeta.naiv)
  
  ##############################covariance matrix of parameter estimation by a given working model
  V.msR <- model.geeglm$geese$vbeta  #robust covariance matrix
  p=dim(V.msR)[1]
  
  ################################model based covairance matrix
  V.modelbase=model.geeglm$geese$vbeta.naiv
  #################################################CIC.H 
  CIC.H <- sum(diag(OmegaI%*%V.msR))
  
  Q=V.msR%*%solve(V.modelbase)
  
  CIC.R=sum(diag(Q))
  ###################################################### 
  C1=sum(diag(Q))/dim(Q)[1]
  C2=sum(diag(Q%*%Q))/dim(Q)[1]
  #estimated mean and observed values
  #RJ=sqrt((1-C1)^2+(1-C2)^2)
  # Delta1=C2-2*C1+1
  # Delta2=sum((log(eigen(Q)$value))^2) 
  mu.R <- model.geeglm$fitted.values
  id=as.numeric(model.geeglm$id)
  m=length(unique(id))
  ni=as.numeric(table(id))
  y <- model.geeglm$y
  X=as.matrix(model.matrix(model.geeglm))
  if(datatype=="continuous"){
    
    phi <-1
    sigma=summary(model.geeglm)$dispersion$Estimate
    
    quasi.R <- sum(-(y-mu.R)^2/(2*phi))   #quasilikelihood for normal model
    
    EQICD=sum((y-mu.R)^2/phi+log(2*pi*sigma))
    
    
    OmegaRI=matrix(0,p,p)
    for(i in 1:m){
      s=ni[i] 
      Di <- matrix(X[id==i,],nr=s) 
      Ai=diag(sigma,nr=s,nc=s)
      OmegaRI= OmegaRI+t(Di)%*%solve(Ai)%*%Di/phi
      
    } 
    
    
  }else if(datatype=="binary"){
    
    phi <- 1  #scale for binary data
    
    quasi.R <- sum(y*log(mu.R/(1-mu.R))+log(1-mu.R)) #quasilikelihood for binomial model
    
    
    EQICD=-2*sum(y*log(mu.R/(1-mu.R))+log(1-mu.R))/phi+sum(log(2*pi*phi*mu.R*(1-mu.R)))
    
    OmegaRI=matrix(0,p,p)
    for(i in 1:m){
      s=ni[i] 
      mu.i=mu.R[id==i]
      Di <- matrix(rep(mu.i*(1-mu.i),p),s,p)*X[id==i,]  
      Ai=diag(mu.i*(1-mu.i),nr=s,nc=s)
      OmegaRI= OmegaRI+t(Di)%*%solve(Ai)%*%Di/phi
      
    } 
    
  }else if(datatype=="count"){
    
    phi =summary(model.geeglm)$dispersion$Estimate
    
    
    quasi.R <- sum(y*log(mu.R)-mu.R)/phi   #quasilikelihood for loglinear model#QIC函数中的拟似然函数没有除以超散布参数phi
    
    #A=diag(mu.R*phi)
    EQICD=sum(y*log(y/mu.R)-(y-mu.R))/phi+sum(log(2*pi*phi*mu.R))
    OmegaRI=matrix(0,p,p)
    for(i in 1:m){
      s=ni[i] 
      mu.i=mu.R[id==i]
      Di <- matrix(rep(mu.i,p),s,p)*X[id==i,]  
      Ai=diag(mu.i,nr=s,nc=s)
      OmegaRI= OmegaRI+t(Di)%*%solve(Ai)%*%Di/phi
      
    } 
    
  }
  
  
  
  
  CIC<- sum(diag(OmegaRI%*%V.msR))
  
  QIC<- (-2)*quasi.R + 2*CIC
  QIC.R <- (-2)*quasi.R + 2*CIC.R   
  QICu=(-2)*quasi.R + 2*p            
  QIC.H <- (-2)*quasi.R + 2*CIC.H    
  EQIC= EQICD+2*CIC
  EQIC.H= EQICD+2*CIC.H
  EQIC.R= EQICD+2*CIC.R
 
  output <- c(QIC,QIC.R,QIC.H,QICu,CIC,CIC.H,CIC.R,EQIC,EQIC.H,EQIC.R)
  
  names(output) <- c('QIC','QIC.R','QIC.H','QICu','CIC',"CIC.H",'CIC.R',"EQIC",'EQIC.H','EQIC.R')
  
  output
  
}


############################################
#Gaussian likelihood criteria
########################################

Gau.geeglm=function(datatype,model.geeglm)
{
  
  mu.R <- as.numeric(model.geeglm$fitted.values)
  
  y <- model.geeglm$y
  
  corr=model.geeglm$corstr
  ni=as.numeric(table(model.geeglm$id))
  m=length(ni)
  M=length(y)
  id=as.numeric(model.geeglm$id)
  alpha=as.numeric(model.geeglm$geese$alpha)
  
  if(datatype=="continuous"){  
    phi <- 1
    sigma=summary(model.geeglm)$dispersion$Estimate
    A=diag(sigma,M)
  }else if(datatype=="binary"){ 
    phi <- 1 #scale for binary data
    A=diag(mu.R*(1-mu.R))
  }else if(datatype=="count"){
    phi <- summary(model.geeglm)$dispersion$Estimate
  #
    A=diag(mu.R)
  }
  
  
  err=y-mu.R
  gau=0
  
  if(corr=="independence"){
    for(i in 1:m){
      s=ni[i]
      erri=err[id==i]
      Ai=matrix(A[id==i,id==i],s,s)
      Vi=phi*Ai
      gau=gau-1/2*(erri%*%solve(Vi)%*%as.matrix(erri,nr=s)+ log(det(2*pi*Vi)) )
      
    }
    
  }else if(corr=="exchangeable"){
    
    for(i in 1:m){
      s=ni[i]
      erri=err[id==i]
      Ai=matrix(A[id==i,id==i],s,s)
      Rex=cs(alpha,s)
      Ai.half=Ai^(1/2)
      Vi=phi*Ai.half%*%Rex%*%Ai.half
      gau=gau-1/2*(erri%*%solve(Vi)%*%as.matrix(erri,nr=s)+ log(det(2*pi*Vi)) )
      
    }
    
  }else if(corr=="ar1"){
    
    for(i in 1:m){
      s=ni[i]
      erri=err[id==i]
      Ai=matrix(A[id==i,id==i],s,s)
      Rar=ar1(alpha,s)
      Ai.half=Ai^(1/2)
      Vi=phi*Ai.half%*%Rar%*%Ai.half
      gau=gau-1/2*(erri%*%solve(Vi)%*%as.matrix(erri,nr=s)+ log(det(2*pi*Vi)) )
      
    }
    
  }
  
  gausslogLK=-2*gau
  return(gausslogLK)
}
#######################Chen and Lazar (2012)  ELR with a stationary structure  ######################


elr.stall <-  function(datatype,model.geeglm) { 
  
  library(emplik)
  
  mu.R <-as.numeric( model.geeglm$fitted.values)
  
  y <- model.geeglm$y
  X=model.geeglm$geese$X
  p=length(model.geeglm$geese$beta)
  corr=model.geeglm$corstr
  id=as.numeric(model.geeglm$id)
  ni=as.numeric(table(id))
  M=length(y)
  nmax=max(ni)
  
  alpha=as.numeric(model.geeglm$geese$alpha)
  
  if(datatype=="continuous"){  
   
    sigma=summary(model.geeglm)$dispersion$Estimate
    A=diag(sigma,M)
    D=array(NA,c(nmax,p,m))
    for(i in 1:m){
      s=ni[i] 
      D[1:s,,i] <- X[id==i,]  
    }
    
  }else if(datatype=="binary"){ 
    
    sigma=mu.R*(1-mu.R)
    A=diag(sigma)
    D=array(NA,c(nmax,p,m))
    for(i in 1:m){
      s=ni[i] 
      mui=mu.R[id==i]
      D[1:s,,i] <-matrix(rep(mui*(1-mui),p),s,p)*X[id==i,]  
    } 
    
  }else if(datatype=="count"){
    sigma=mu.R
    A=diag(sigma)
    D=array(NA,c(nmax,p,m))
    for(i in 1:m){
      s=ni[i] 
      D[1:s,,i] <-matrix(rep(mu.R[id==i],p),s,p)*X[id==i,]  
    } 
  } 
  
  
  
  error=y-mu.R
  perror=error/sqrt(sigma)
  phi <- sum(perror^2)/(M-p)
  
  
  r <- p+nmax-1   
  g <- matrix(0, r,m)
  
  for (i in 1:m){
    s <- ni[i]
    if(corr=="independence"){
      alpha.EL=rep(0,s-1)
    }
    else if(corr=="exchangeable"){
      alpha.EL=rep(alpha,s-1)      
    }else if(corr=="ar1"){
      alpha.EL=rep(alpha,s-1)^c(1:(s-1))     
    }
    
    Di <- matrix(D[1:s,,i],nr=s,nc=p)
    Ri <- top(alpha.EL)
    error.i=error[id==i]
    Ai=A[id==i,id==i]
    Ai.half=Ai^(1/2)
    Vi=phi*Ai.half%*%Ri%*%Ai.half
    p.erri <- perror[id==i]
    
    g[1:p,i] <- t(Di)%*%solve(Vi)%*%matrix(error.i,nr=s)  ## the first p components of g
    if(s==1){g[-c(1:p),i]=0}
    else if(s!=1){
      for (k in 1:(s-1)){
        sum <- 0
        for (j in 1:(s-k)){
          sum <- sum + p.erri[j]*p.erri[j+k]
        }
        g[k+p,i] <- sum-alpha.EL[k]*(s-k-p/m)*phi
      }     
    }  
  }
  g <- t(g)
  g.mu <- rep(0,r)
  
  el.test(g, g.mu,gradtol=1e-9)$"-2LLR"
}





elr.st <-  function(xmat,ymat,id,beta,alpha,phi,datatype) 
{  
  ni=as.numeric(table(id)) # the number of observations for each subject
  m=length(unique(id))  #sample size
  nmax=max(ni)
  M=sum(ni)
  p=length(beta)
  g <- matrix(0, p+nmax-1,m)
  R <- top(alpha)
  
  if(datatype=="continuous")
  {
    for (i in 1:m) {
      s=ni[i]
      fitted= xmat[1:s,,i] %*% beta
      D<-matrix(xmat[1:s,,i],nc=p)
      A<-diag(1,s)
      error <- ymat[1:s,i]-fitted
      p.err<-  error
      
      RR <- matrix(R[1:s,1:s],nc=s)
      A.half <- A^(1/2)
      g[1:p,i] <- crossprod(D, A.half %*% solve(RR, A.half %*% error))*phi^(-1)# the first p components of g (gee)
      if(ni[i]==1){stop} else{
        for(k in 1:(s-1)){
          sum=0
          for(j in 1:(s-k)){
            sum=sum+p.err[j]*p.err[j+k]
            
          }
          g[p+k,i]=sum-alpha[k]*(s-k-p/m)*phi 
          
        }
      }
    }
    
  }else if(datatype=="binary"){
    for (i in 1:m) {
      s=ni[i]
      fitted <-  plogis( xmat[1:s,,i] %*% beta )
      error <- ymat[1:s,i]-fitted
      p.err<- error*(fitted*(1-fitted))^(-1/2) 
      D <- matrix(rep(fitted*(1-fitted),p),s,p)*xmat[1:s,,i]
      A <- diag(as.vector((fitted*(1-fitted))^(-1)),nc=s)
      
      error <- ymat[1:s,i]-fitted
      
      
      RR <- matrix(R[1:s,1:s],nc=s)
      A.half <- A^(1/2)
      g[1:p,i] <- crossprod(D, A.half %*% solve(RR, A.half %*% error))*phi^(-1)# the first p components of g (gee)
      if(ni[i]==1){stop} else{
        for(k in 1:(s-1)){
          sum=0
          for(j in 1:(s-k)){
            sum=sum+p.err[j]*p.err[j+k]
            
          }
          g[p+k,i]=sum-alpha[k]*(s-k-p/m)*phi 
          
        }
      }
    }
    
    
    
  } else if (datatype=="count")
  {
    for (i in 1:m) 
    {
      s=ni[i]
      fitted <- exp( xmat[1:s,,i] %*% beta)
      error <- ymat[1:s,i]-fitted
      p.err <- error*fitted^(-1/2)
      D <- matrix(rep(c(fitted),p),s,p)*xmat[1:s,,i]
    
      A <- diag(as.vector(fitted^(-1)),nc=s)
      error <- ymat[1:s,i]-fitted
      RR <- matrix(R[1:s,1:s],nc=s)
      A.half <- A^(1/2)
      g[1:p,i] <- crossprod(D, A.half %*% solve(RR, A.half %*% error))*phi^(-1)# the first p components of g (gee)
      if(ni[i]==1){stop} else{
        for(k in 1:(s-1)){
          sum=0
          for(j in 1:(s-k)){
            sum=sum+p.err[j]*p.err[j+k]
            
          }
          g[p+k,i]=sum-alpha[k]*(s-k-p/m)*phi 
          
        }
      }
    }
  }
  g <- t(g)
  g.mu <- rep(0,p+nmax-1)
  el.test(g, g.mu,gradtol=1e-9)$"-2LLR"
}
#####################################################################
########################################################
#empirical likelihood method
#########################################################

#################################
# Gosho et al 2011 method
##################################
shof=function( datatype,model.geeglm)
{
  require(magic)
  require(Rlab)
  y <- model.geeglm$y
  x=as.matrix(model.matrix(model.geeglm))
  mu.R <-as.numeric( model.geeglm$fitted.values)
  corr=model.geeglm$corstr
  ni=as.numeric(table(model.geeglm$id))
  n=unique(ni)
  m=length(ni)
  p=dim(x)[2]
  M=length(y)
  id=as.numeric(model.geeglm$id)
  a=as.numeric(model.geeglm$geese$alpha)
  phi.i=summary(model.geeglm)$dispersion$Estimate
  
  if (corr=="independence") {R <- diag(n)
  }  else if (corr=="exchangeable"){
    R <- cs(a,n)
  }  else if (corr=="ar1"){
    R <- ar1(a,n)
  } else if (corr=="ma1"){
    R <- ma1(a,n)
  }  else if (corr=="stat"){
    R <- top(a)
  }
  X <- aperm(array(t(x),c(p,n,m)),c(2,1,3))           ## data transformation
  Y <-  matrix(y, n,m)
  
  crs=matrix(0,n,n)
  vr=matrix(0,n,n)
  for (i in 1:m){
   
    if(datatype=="continuous"){
     
      D<-matrix(X[,,i],nc=p)
      A<-diag(1,n)
      error <- Y[,i]-mu.R[id==i]
      
    }else if(datatype=="count"){
  
     fitted <- mu.R[id==i]
    error <- Y[,i]-fitted
    D <- matrix(rep(fitted,p),n,p)*X[,,i]
    A <- diag(as.vector((fitted)))
    } else if (datatype=="binary"){
      fitted <-  plogis( mu.R[id==i] )
      error <- Y[,i]-fitted
       
      D <- matrix(rep(fitted*(1-fitted),p),s,nc=p)*X[,,i]
      A <- diag(as.vector((fitted*(1-fitted))^(-1)))  
      
    }
    
    Ai.half <- A^(1/2)
    
    vi=phi.i*Ai.half%*%R%*%Ai.half
    
    crs=crs+error%*%t(error)
    vr=vr+vi
  }
  
  AA=crs%*%solve(vr)-diag(n)
  CR=sum(diag(AA%*%AA))
  
  return(CR)
}
