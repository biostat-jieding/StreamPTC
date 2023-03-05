

#==========================================================================#
# Promotion time cure model (PTC) -- based on: NPMLE from Portier et al. 2017
#==========================================================================#

#==== Function for fitting PTC with NPMLE ====#
#' @title Promotion time cure model based on the nonparametric maximum likelihood estimation procedure
#'
#' @description Fit promotion time cure model based on the nonparametric maximum likelihood estimation (NPMLE) procedure.
#' More details about the classical NPMLE method can be found in Portier et al. (2017). 
#' 
#' @aliases PTC.NPMLE.fit
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' @param boots a logical value. The default specification is \code{TRUE}, indicating that the standard errors will be provided using bootstrap.
#' @param nboot specifies the number of bootstrap sampling. The default \code{nboot = 100}.
#' @param trace a binary variable that determines whether to display the information for each bootstrap iteration.
#' @param maxit specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{maxit = 5e3}. 
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param simplify a logical value. The default specification is \code{TRUE}, indicating that whether we will use the internal \code{optim} function in \code{R} to solve the method.
#' 
#' @return The fitted results are returned (a list). The estimates of regression coefficients is \code{res}.
#' 
#' @references PORTIER, F., EL GHOUCH, A. and VAN KEILEGOM, I. (2017). Efficiency and bootstrap in the promotion time cure model. Bernoulli 23 3437–3468.
#' 
#' @examples 
#' ### ==== An example for fitting the promotion time cure model ==== ###
#' 
#' ## ---- generate the simulated dataset ---- ##
#' sdata.full  <- sdata.PTC(N=200)
#' yobs  <- sdata.full$sdata$yobs
#' delta <- sdata.full$sdata$delta
#' X     <- as.matrix(sdata.full$sdata$X)
#' rm(sdata.full)
#' 
#' ## ---- fit the promotion time cure model ---- ##
#' # estimator of the regression coefficients
#' fit.PTC <- PTC.NPMLE.fit(yobs,delta,X)
#' print(fit.PTC$res)
#' # estimator of the baseline distribution at certain points
#' tm <- c(0.1,0.3,0.5,0.7,0.9)
#' Est.Lam0t <- PTC.NPMLE.Lam(tm=tm,bet=fit.PTC$res[,1],yobs=yobs,delta=delta,X=X)
#' print(data.frame(Est=Est.Lam0t,row.names=tm))
#' 
#' @export PTC.NPMLE.fit
PTC.NPMLE.fit <-  function(yobs,delta,X,
                           LinkFunc=list(D0=exp,D1=exp),
                           boots=FALSE,nboot=100,
                           eps=1e-6,maxit=5e3,
                           simplify=TRUE){
  # need: survival
  
  ### Preparations
  N <- length(yobs)
  XI <- cbind(1,X)
  pbet <- ncol(XI)

  ### calculate initial values for bet: by Cox model (bet[1] by cure rate from KM estimation at the last event time)
  KM.fit <- survival::survfit(survival::Surv(yobs,delta)~1)
  cure.rate <- min(KM.fit$surv)  
  bet.init <- c(ifelse(cure.rate==0,0,log(-log(cure.rate))),rep(0,pbet-1))
  
  ### calculate MLE for beta
  if(simplify==TRUE){
    
    numit <- 1
    bet.old <- bet.init
    repeat{
      InfoMScore <- PTC.NPMLE.ScoreInfoM(bet=bet.old,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=TRUE)
      Score <- InfoMScore$Score
      InfoM <- InfoMScore$InfoM
      dev <- MASS::ginv(InfoM)%*%Score # solve(InfoM,Score)
      bet <- bet.old + as.vector(dev)
      if( max(abs(bet-bet.old))>eps & numit<maxit ){
        bet.old <- bet
        numit <- numit + 1
      }else{
        break
      }
    }
    convergence <- numit<maxit
    
  }else{
    
    sol <- stats::optim(par=bet.init,fn=PTC.NPMLE.Beta.LogLik,
                        control=list(maxit=maxit,fnscale=-1,reltol=eps),
                        yobs=yobs,delta=delta,XI=XI) # sol
    bet <- sol$par
    convergence <- sol$convergence
    
  }
  
  
  ### calculate SEs
  if(boots==FALSE){
    
    ## calculate SEs for bet using explicit formula !!!
    InfoMScore <- PTC.NPMLE.ScoreInfoM(bet=bet,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=TRUE)
    InfoM <- InfoMScore$InfoM
    Score <- InfoMScore$Score
    V <- solve(InfoM)
    bet.se <- sqrt(diag(V)/N)
    bet.boots <- NULL

  }else if(boots==TRUE){
    
    ## calculate SEs for bet using bootstrap (direct or weighted)
    bet.boots <- array(0,dim=c(nboot,pbet))
    for(iboot in 1:nboot){
      # bootstrap weights
      e <- rexp(N)
      w <- e/mean(e)
      # solve the current estimator
      sol.iboot <- stats::optim(par=bet,fn=PTC.NPMLE.Beta.LogLik,
                                control=list(maxit=1e4,fnscale=-1),
                                yobs=yobs,delta=delta,XI=XI,w=w)
      bet.boots[iboot,] <- sol.iboot$par
    }
    bet.boots.center <- bet.boots - matrix(apply(bet.boots,2,mean),nrow=nboot,ncol=pbet,byrow=TRUE)
    V.bet <- (t(bet.boots.center) %*% bet.boots.center)/nboot
    bet.se <- sqrt(diag(V.bet))
    V <- V.bet * N
    InfoM <- MASS::ginv(V,tol=1e-8)
  }

  
  ### tidy the results: inference
  zvalue.bet <- bet/bet.se
  pvalue.bet <- 2*(1-pnorm(abs(zvalue.bet)))
  res <- data.frame(Est=bet,SE=bet.se,zvalue=zvalue.bet,pvalue=pvalue.bet,
                    row.names=c("(Intercept)",colnames(X)))
  
  ### output
  out <- list(
    convergence = convergence,
    res=res,
    InfoM=InfoM,
    Score=Score,
    bet.boots=bet.boots
  )
  return(out)
  
}


#==== Loss function for bet (Profiled log-likelihood function) ====#
#' @title Loss function for bet in NPMLE
#'
#' @description  Loss function for the promotion cure model (Profiled log-likelihood function in NPMLE).
#' 
#' @aliases PTC.NPMLE.Beta.LogLik
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param XI a matrix of covariates with 1 as the first column, that is, \code{XI=cbind(1,X)}.
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' @param w weights for the weighted bootstrap method. The default is \code{NULL}, and it will be treatead as equal weights (\code{w=rep(1,length(yobs))}).
#' 
#' @export PTC.NPMLE.Beta.LogLik
PTC.NPMLE.Beta.LogLik <- function(bet,yobs,delta,XI,LinkFunc=list(D0=exp,D1=exp),w=NULL){
  # for maximization
  
  # prepare: etaXIbet and Rbet
  N <- length(yobs)
  if(is.null(w)){w <- rep(1,N)}
  etaXIbet <- as.vector(LinkFunc$D0(XI %*% bet))
  dRbet <- sapply(1:N,function(i){
    delta[i]*sum( etaXIbet*(yobs>=yobs[i])*w )/N
  }) # = Rbet * delta
  # Delta <- 1*(yobs<=max(yobs[delta==1]))
  # Rbet <- sapply(yobs,function(yobsi){
  #   sum( etaXIbet*(Delta*(yobs>=yobsi)+1-Delta)*w )/N
  # })
  
  # calculate lambda
  Rbet.min <- min(dRbet[delta==1])
  interval <- c(Rbet.min-sum(delta*w)/N,Rbet.min-min(w[delta==1])/N)
  lambet <- uniroot(PTC.NPMLE.Lambda.Solve,interval,tol = .Machine$double.eps^0.75,
                    dRbet=dRbet,delta=delta,w=w)$root
  
  # calculate the loss function for beta (log form)
  val.log <- sum( log(etaXIbet[delta==1]/(dRbet[delta==1]-lambet))*w[delta==1] ) - N*lambet
  
  # output
  return(val.log)
  
}


#==== Nonparametric component ====#
#' @title Nonparametric component in NPMLE
#'
#' @description  Nonparametric component for the promotion cure model at specified time points.
#' 
#' @aliases PTC.NPMLE.Lam
#'
#' @param tm The time points that the nonparametric baseline distribution will be estimated at. The default specification is \code{TRUE}.
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param XI a matrix of covariates with 1 as the first column, that is, \code{XI=cbind(1,X)}.
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' @param w weights for the weighted bootstrap method. The default is \code{NULL}, and it will be treatead as equal weights (\code{w=rep(1,length(yobs))}).
#' 
#' @export PTC.NPMLE.Lam
PTC.NPMLE.Lam <- function(tm,bet,yobs,delta,X,LinkFunc=list(D0=exp,D1=exp),w=NULL){
  
  # prepare: etaXIbet and Rbet
  N <- length(yobs)
  XI <- cbind(1,X)
  if(is.null(w)){w <- rep(1,N)}
  etaXIbet <- as.vector(LinkFunc$D0(XI %*% bet))
  dRbet <- sapply(1:N,function(i){
    delta[i]*sum( etaXIbet*(yobs>=yobs[i])*w )/N
  }) # = Rbet * delta
  
  # calculate lambda
  Rbet.min <- min(dRbet[delta==1])
  interval <- c(Rbet.min-sum(delta*w)/N,Rbet.min-min(w[delta==1])/N)
  lambet <- uniroot(PTC.NPMLE.Lambda.Solve,interval,tol = .Machine$double.eps^0.75,
                    dRbet=dRbet,delta=delta,w=w)$root
  
  # calculate the Lam(t) at specified time points tm
  Lam <- sapply(tm,function(tmj){
    sum( (w*(yobs<=tmj))[delta==1]/(dRbet[delta==1]-lambet) ) / N
  })
  
  # output
  return(Lam)
  
}


#==== An auxiliary function for solving lambda given bet ====#
#' @title An auxiliary function for solving lambda given bet
#'
#' @description  An auxiliary function for solving lambda given bet.
#' 
#' @aliases PTC.NPMLE.Lambda.Solve
#'
#' @param lam unknown parameter corresponding to the largrane multiplier.
#' @param Rbet a matrix needed in this function.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param w weights for the weighted bootstrap method. The default is \code{NULL}, and it will be treatead as equal weights (\code{w=rep(1,length(yobs))}).
#' 
#' @export PTC.NPMLE.Lambda.Solve
PTC.NPMLE.Lambda.Solve <- function(lam,dRbet,delta,w=NULL){
  
  if(is.null(w)){w <- rep(1,length(delta))}
  val <- sum(w[delta==1]/(dRbet[delta==1]-lam))/length(delta) - 1
  return(val)
  
}

#==== Obtain score vector and information matrix ====#
#' @title Score vector and information matrix in NPMLE
#'
#' @description calculate score vector and information matrix in NPMLE.
#' 
#' @aliases PTC.NPMLE.ScoreInfoM
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' @param IsScore whether the score vector is calculated or not.
#' @param IsInfoM whether the information matrix is calculated or not.
#' 
#' @export PTC.NPMLE.ScoreInfoM
PTC.NPMLE.ScoreInfoM <-  function(bet,yobs,delta,X,LinkFunc=list(D0=exp,D1=exp),
                                  IsScore=FALSE,IsInfoM=TRUE){
  
  # prepare elements
  XI <- cbind(1,X)
  N <- length(yobs)
  yobsd <- yobs[delta==1]
  ord <- order(yobsd)
  XIbet <- as.vector(XI %*% bet)
  etaXIbet <- LinkFunc$D0(XIbet)
  eta1XIbet <- LinkFunc$D1(XIbet)
  tau <- max(yobs[delta==1]) # ifelse(is.null(tau),max(yobs[delta==1]),max(tau,max(yobs[delta==1])))
  Delta <- 1*(yobs<=tau)
  d <- XI*eta1XIbet/etaXIbet
  dd <- d[delta==1,,drop=F]
  dRbet <- sapply(1:N,function(i){
    delta[i]*sum( etaXIbet*(yobs>=yobs[i]) )/N
  }) # = Rbet * delta
  Rbetd <- dRbet[delta==1]
  # calculate lambda
  Rbet.min <- min(dRbet[delta==1])
  interval <- c(Rbet.min-sum(delta)/N,Rbet.min-1/N)
  lambet <- uniroot(PTC.NPMLE.Lambda.Solve,interval,tol = .Machine$double.eps^0.75,
                    dRbet=dRbet,delta=delta)$root
  # some other elements
  Dd <- t(sapply(1:length(yobsd),function(j){apply(d*(etaXIbet*(Delta*(yobs>=yobsd[j])+1-Delta)),2,mean)}))
  ch2 <- sum(1/(Rbetd*(Rbetd-lambet)))/N
  ch1 <- apply(Dd/(Rbetd*(Rbetd-lambet)),2,sum)/N
  ch <- ch1/ch2
  hd <- t(t(Dd)-ch)/Rbetd
  
  out <- list()
  ## prepare score vector
  if(IsScore==TRUE){
    Score1 <- apply(dd-hd,2,sum) / N
    Score2 <- ch 
    Score <- (Score1 - Score2) # == apply(dd-Dd/(Rbetd-lambet),2,sum)/N
    out <- c(out,list(Score=Score))
  }
  ## prepare information matrix
  if(IsInfoM==TRUE){
    ## Way 1: the covariance of score vector (Beyhum et al., 2017)
    # InfoM <- t(dd-hd)%*%(dd-hd)/N
    # Way 2: the covariance of score vector (Portier et al., 2017)
    InfoM <- 0
    for(j in 1:length(yobsd)){
      InfoM <- InfoM +
        t( t(t(d)-hd[j,])*(etaXIbet*(Delta*(yobs>=yobsd[j])+1-Delta)) )%*%t(t(d)-hd[j,]) /
        N^2 / (Rbetd[j]-lambet)
    }
    ## Way 3: the derivative of score vector
    # lambet.D1 <- apply(Dd/(Rbetd-lambet)^2,2,sum)/sum(1/(Rbetd-lambet)^2)
    # eta2XIbet <- LinkFunc$D2(XIbet)
    # eta210 <- (eta2XIbet/etaXIbet-(eta1XIbet/etaXIbet)^2)
    # Ebet1 <- t(XI)%*%(XI*eta210*delta)/N
    # Rbetd.D2 <- lapply(1:length(yobsd),function(i){
    #   t(XI) %*% (XI*eta2XIbet*(Delta*(yobs>=yobsd[i])+1-Delta)) / N
    # })
    # Ebet2d <- lapply(1:length(yobsd),function(i){
    #   Rbetd.D2[[i]]/(Rbetd[i]-lambet) - Dd[i,]%*%t(Dd[i,]-lambet.D1)/(Rbetd[i]-lambet)^2
    # })
    # Ebet2 <- 0
    # for(j in 1:length(yobsd)){
    #   Ebet2 <- Ebet2 + Ebet2d[[j]]/N
    # }
    # InfoM <- Ebet2-Ebet1
    # out
    out <- c(out,list(InfoM=InfoM))
  }
  
  ## output
  return(out)
  
}


#==== Survival function for the whole population (for promotion time cure model) ====#
#' @title Survival function for the whole population (for promotion time cure model)
#'
#' @description Calculate the survival function for the whole population (for promotion time cure model) at various time points.
#' 
#' @aliases PTC.NPMLE.Stx 
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param bet unknown parameters corresponding to latency part,   or the fitted values of regression coefficients in the latency part returned by the function \code{\link{SMC.PH.fit()}}.
#' @param cross a logical value. The default is \code{FALSE}. 
#' @param tm the time points that the survival function will be calculated at. The default is \code{NULL}.
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' 
#' @export PTC.NPMLE.Stx 
PTC.NPMLE.Stx <- function(yobs,delta,X,bet,cross=FALSE,tm=NULL,LinkFunc=list(D0=exp,D1=exp)){ # for uncured 
  
  # calculate
  XI <- cbind(1,X)
  etaXbet <- as.vector(LinkFunc$D0(XI%*%bet))
  if(cross){
    # calculate baseline cumulative function
    Lamt <- PTC.NPMLE.Lam(tm=tm,bet=bet,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc) # Lamt[tm>max(yobs[delta==1])] <- 1; Lamt[tm<min(yobs[delta==1])] <- 0
    # calculate Survivals
    Stx <- exp(-outer(etaXbet,Lamt,function(x,y){x*y}))
  }else{
    # calculate baseline cumulative function
    Lamt <- PTC.NPMLE.Lam(tm=yobs,bet=bet,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc)
    # calculate Survivals
    Stx <- exp(-etaXbet*Lamt)
  }
  # output
  return(list(Stx=Stx))
  
}

#==== Obtain half of the information matrix ====#
#' @title Half of the information matrix in NPMLE
#'
#' @description calculate the half of the information matrix in NPMLE
#' 
#' @aliases PTC.NPMLE.InfoM.Half
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' 
#' @export PTC.NPMLE.InfoM.Half
PTC.NPMLE.InfoM.Half <-  function(bet,yobs,delta,X,LinkFunc=list(D0=exp,D1=exp)){
  
  # prepare elements
  XI <- cbind(1,X)
  N <- length(yobs)
  pbet <- length(bet)
  yobsd <- yobs[delta==1]
  idx.delta <- which(delta==1)
  ord <- order(yobsd)
  XIbet <- as.vector(XI %*% bet)
  etaXIbet <- LinkFunc$D0(XIbet)
  eta1XIbet <- LinkFunc$D1(XIbet)
  tau <- max(yobs[delta==1]) # ifelse(is.null(tau),max(yobs[delta==1]),max(tau,max(yobs[delta==1])))
  Delta <- 1*(yobs<=tau)
  d <- XI*eta1XIbet/etaXIbet
  dd <- d[delta==1,,drop=F]
  dRbet <- sapply(1:N,function(i){
    delta[i]*sum( etaXIbet*(yobs>=yobs[i]) )/N
  }) # = Rbet * delta
  Rbetd <- dRbet[delta==1]
  # calculate lambda
  Rbet.min <- min(dRbet[delta==1])
  interval <- c(Rbet.min-sum(delta)/N,Rbet.min-1/N)
  lambet <- uniroot(PTC.NPMLE.Lambda.Solve,interval,tol = .Machine$double.eps^0.75,
                    dRbet=dRbet,delta=delta)$root
  # some other elements
  Dd <- t(sapply(1:length(yobsd),function(j){apply(d*(etaXIbet*(Delta*(yobs>=yobsd[j])+1-Delta)),2,mean)}))
  ch2 <- sum(1/(Rbetd*(Rbetd-lambet)))/N
  ch1 <- apply(Dd/(Rbetd*(Rbetd-lambet)),2,sum)/N
  ch <- ch1/ch2
  hd <- t(t(Dd)-ch)/Rbetd
  
  ## prepare information matrix (half)
  InfoM.Half <- array(0,dim=c(N^2,pbet))
  for(j in 1:length(yobsd)){
    InfoM.Half[(idx.delta[j]-1)*N+(1:N),] <- 
      t(t(d)-hd[j,]) * sqrt((etaXIbet*(Delta*(yobs>=yobsd[j])+1-Delta))/(Rbetd[j]-lambet))/N
  }
  
  ## output
  return(InfoM.Half)
  
}


#==== Obtain individual level score vecto ====#
#' @title Individual level score vector in NPMLE
#'
#' @description calculate the individual level score vector in NPMLE
#' 
#' @aliases PTC.NPMLE.Score.Individual
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' @param Iscov whether convert the individual level form into a covariance form.
#' 
#' @export PTC.NPMLE.Score.Individual
PTC.NPMLE.Score.Individual <-  function(bet,yobs,delta,X,LinkFunc=list(D0=exp,D1=exp),Iscov=TRUE){
  # Score is the first  derivative of [positive] log-likelihood
  
  ## prepare elements
  N <- length(yobs)
  XI <- cbind(1,X)
  pbet <- ncol(XI)
  yobsd <- yobs[delta==1]
  ord <- order(yobsd)
  XIbet <- as.vector(XI %*% bet)
  etaXIbet <- LinkFunc$D0(XIbet)
  eta1XIbet <- LinkFunc$D1(XIbet)
  tau <- max(yobs[delta==1]) # ifelse(is.null(tau),max(yobs[delta==1]),max(tau,max(yobs[delta==1])))
  Delta <- 1*(yobs<=tau)
  d <- XI*eta1XIbet/etaXIbet
  dd <- d[delta==1,,drop=F]
  dRbet <- sapply(1:N,function(i){
    delta[i]*sum( etaXIbet*(yobs>=yobs[i]) )/N
  }) # = Rbet * delta
  Rbetd <- dRbet[delta==1]
  # calculate lambda
  Rbet.min <- min(dRbet[delta==1])
  interval <- c(Rbet.min-sum(delta)/N,Rbet.min-1/N)
  lambet <- uniroot(PTC.NPMLE.Lambda.Solve,interval,tol = .Machine$double.eps^0.75,
                    dRbet=dRbet,delta=delta)$root
  # some other elements
  Dd <- t(sapply(1:length(yobsd),function(j){apply(d*(etaXIbet*(Delta*(yobs>=yobsd[j])+1-Delta)),2,mean)}))
  ch2 <- sum(1/(Rbetd*(Rbetd-lambet)))/N
  ch1 <- apply(Dd/(Rbetd*(Rbetd-lambet)),2,sum)/N
  ch <- ch1/ch2
  hd <- t(t(Dd)-ch)/Rbetd
  # calculate U
  Ud <- t(t(dd-hd)-ch) # here "minus ch" is wrong, we should also expand ch
  U <- array(0,dim=c(N,pbet))
  U[delta==1,] <- Ud
  if(Iscov==TRUE){
    out <- t(U)%*%U/N
  }else{
    out <- U
  }
  
  ## output
  return(out)
  
}

