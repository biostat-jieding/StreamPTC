#------------------------------------------------------------------------#
# data generating function (low-dimension) ####
#------------------------------------------------------------------------#
#' @title The data generating function under the promotion time cure model (low-dimensional)
#'
#' @description Generate a simulated dataset from the promotion time cure model with low dimensional covariates.
#' 
#' @aliases sdata.PTC
#'
#' @param N the sample size.
#' @param bet true value of coefficients (include the intercept). Certain default value is set.
#' @param Lam0t true baseline distribution function. Certain default value is set.
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' @param avalue a value used to control the cure rate.
#' @param cvalue a value used to control the censoring rate.
#' @param cureRate desired cure rate. 
#' @param censoringRate desired censoring rate.
#' @param Interaction whether to include the interaction term (then bet should be of 4-dimension).
#'   
#' @export sdata.PTC
sdata.PTC <- function(N,
                      bet = c(-1,1,1),
                      Lam0t=function(t){ifelse(t>1,1,ifelse(t<0,0,t))},
                      LinkFunc=list(D0=exp,D1=exp),
                      avalue=0,cvalue=2.5,cureRate=NULL,censoringRate=NULL,
                      Interaction=FALSE
){
  # Generate Data from Ptomotion Time Cure Model: S(t|x)=exp(-eta(Xbet)Lam0t)

  ### basic setup for several elements and functions
  Ftx <- function(t,x,bet,Lam0t,LinkFunc){ # pupulation distribution function
    FF <- 1 - exp( -LinkFunc$D0(sum(c(1,x)*bet))*Lam0t(t) )
    return(FF)
  }
  Covariate <- function(N,avalue,Interaction){ # generate covariates
    X <- array(NA, dim=c(N,2))
    X[,1] <- runif(N,avalue,1+avalue) 
    X[,2] <- rbinom(N,1,0.5)
    if(Interaction==TRUE){
      X <- cbind(X,X[,1]*X[,2])
    }
    return(X)
  }   
  
  ### generate covariates and cure status for patients
  # generate covariate
  X <- Covariate(N,avalue,Interaction)
  # generate cure status
  p.cure <- as.vector(exp(-LinkFunc$D0(cbind(1,X)%*%bet)))
  repeat{ 
    cure.state <- rbinom(N, 1, p.cure)   # 1=cured or 0=uncured
    cure.rate <- mean(cure.state)
    if( is.null(cureRate) ){ break } 
    if( cure.rate <= cureRate + 0.04 && 
        cure.rate >= cureRate - 0.04 ){ break }
  }
  
  # generate survival time
  get.stime <- function(t,u,x,bet,Lam0t,LinkFunc){ u - Ftx(t,x,bet,Lam0t,LinkFunc) }
  stime <- rep(NA, N)
  for(i in 1:N){
    if(cure.state[i]==0){
      stime[i] <- uniroot(get.stime,c(0,1),u=runif(1,0,1-p.cure[i]),x=X[i,],
                          bet=bet,Lam0t=Lam0t,LinkFunc=LinkFunc)$root
    }else{
      stime[i] <- Inf
    }
  }
  
  ### generate censoring time
  repeat{ # censoring.rate is random, to control it 
    ctime <- rexp(N,1/cvalue) # runif(N,0,cvalue) 
    delta <- as.numeric(stime<=ctime)
    censoring.rate <- 1 - mean(delta)  
    if( is.null(censoringRate) ){ break } 
    if( censoring.rate <= censoringRate + 0.05 && 
        censoring.rate >= censoringRate - 0.05 ){ break } 
  }
  
  # delta and observed failure time
  yobs <- pmin(stime,ctime)
  
  # info
  info <- list(
    cure.rate = mean(cure.state),
    censoring.rate = 1-mean(delta)
  ); info
  
  # X <- as.matrix(data.frame(X))
  
  # output
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=data.frame(X)),
    info=info
  )
  return(out)
  
}

#------------------------------------------------------------------------#
# data generating function (high-dimension) ####
#------------------------------------------------------------------------#
#' @title The data generating function under the promotion time cure model (high-dimensional)
#'
#' @description Generate a simulated dataset from the promotion time cure model with relatively high dimensional covariates.
#' 
#' @aliases sdata.PTC.High
#'
#' @param N the sample size.
#' @param bet.noInt true value of coefficients (not include the intercept). Certain default value is set.
#' @param pX the total number of generated covariates. The default value is 10.
#' @param X.type the type of covariates. The default type is "continuous" (normal).
#' @param avalue a value used to control the cure rate. It is also the intercept term.
#' @param cvalue a value used to control the censoring rate.
#' @param cureRate desired cure rate. 
#' @param censoringRate desired censoring rate.
#' @param corstr type of correlation structure.
#' @param rho the degree of correlation.
#' @param tau the truncated point.
#' @param Lam0t true baseline distribution function. Certain default value is set.
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#'   
#' @export sdata.PTC.High
sdata.PTC.High <- function(N,
                           bet.noInt=c(0.6,0.5,0.4),
                           pX=10,
                           X.type="continuous",
                           avalue=-0.5,cvalue=2.2,cureRate=NULL,censoringRate=NULL,
                           corstr=c("IND","EX","AR1")[3],rho=0.4,tau=2,
                           Lam0t=function(t){ifelse(t>1,1,ifelse(t<0,0,t))},
                           LinkFunc=list(D0=exp,D1=exp)
){
  # Generate Data from Ptomotion Time Cure Model: S(t|x)=exp(-eta(Xbet)Lam0t)

  ### basic setup for several elements and functions
  Ftx <- function(t,x,bet,Lam0t,LinkFunc){ # pupulation distribution function
    FF <- 1 - exp( -LinkFunc$D0(sum(c(1,x)*bet))*Lam0t(t) )
    return(FF)
  }
  
  ### generate covariates and cure status for patients
  # generate covariate
  bet <- c(avalue,bet.noInt,rep(0,pX-length(bet.noInt)))
  if(X.type=="continuous"){
    X <- array(NA, dim=c(N,pX))
    if(corstr=="EX"){
      Sigma <- rho*(1-diag(pX))+diag(pX)
    }else if(corstr=="AR1"){
      Sigma <- rho^(abs(row(diag(pX))-col(diag(pX))))
    }else if(corstr=="IND"){
      Sigma <- diag(pX)
    }
    X[,1:pX] <- mvtnorm::rmvnorm(N,mean=rep(0,nrow(Sigma)),sigma=Sigma)
    X <- apply(X,2,function(x){ x[x>tau] <- tau;x[x< -tau]<- -tau; x })
  }else if(X.type=="binary"){
    X <- array(rbinom(N*pX,1,0.5),dim=c(N,pX))
  }
  # generate cure status
  p.cure <- as.vector(exp(-LinkFunc$D0(cbind(1,X)%*%bet)))
  repeat{ 
    cure.state <- rbinom(N, 1, p.cure)   # 1=cured or 0=uncured
    cure.rate <- mean(cure.state)
    if( is.null(cureRate) ){ break } 
    if( cure.rate <= cureRate + 0.04 && 
        cure.rate >= cureRate - 0.04 ){ break }
  }
  
  # generate survival time
  get.stime <- function(t,u,x,bet,Lam0t,LinkFunc){ u - Ftx(t,x,bet,Lam0t,LinkFunc) }
  stime <- rep(NA, N)
  for(i in 1:N){
    if(cure.state[i]==0){
      stime[i] <- uniroot(get.stime,c(0,1),u=runif(1,0,1-p.cure[i]),x=X[i,],
                          bet=bet,Lam0t=Lam0t,LinkFunc=LinkFunc)$root
    }else{
      stime[i] <- Inf
    }
  }
  
  ### generate censoring time
  repeat{ # censoring.rate is random, to control it 
    ctime <- rexp(N,1/cvalue) # runif(N,0,cvalue) 
    delta <- as.numeric(stime<=ctime)
    censoring.rate <- 1 - mean(delta)  
    if( is.null(censoringRate) ){ break } 
    if( censoring.rate <= censoringRate + 0.05 && 
        censoring.rate >= censoringRate - 0.05 ){ break } 
  }
  
  # delta and observed failure time
  yobs <- pmin(stime,ctime)
  
  # info
  info <- list(
    cure.rate = mean(cure.state),
    censoring.rate = 1-mean(delta)
  ); info
  
  # X <- as.matrix(data.frame(X))
  
  # output
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=data.frame(X)),
    info=info
  )
  return(out)
  
}


