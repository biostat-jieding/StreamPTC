
#==========================================================================#
# High-dimensional Variable Selection in Promotion time cure model (ALasso)
#==========================================================================#

#==== High-dimensional Variable Selection in PTC with NPMLE (ALasso, Offline) ====#
#' @title Variable selection in promotion time cure model (via adaptive Lasso)
#' based on the classical nonparametric maximum likelihood estimator and adaptive Lasso penalty
#'
#' @description Conduct variable selection using promotion time cure model based on the nonparametric maximum likelihood estimation (NPMLE) procedure and adaptive Lasso penalty.
#' In other words, this function can analyze survival data with a cure fraction and selection the underlying useful variables.
#'
#' @aliases PTC.NPMLE.ALasso.fit
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' @param tunpara.num the number of tuning parameter values that will be fitted at.
#' @param tunpara.min.ratio Smallest value for tuning parameter, as a fraction of tunpara.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero).
#' @param tunpara A user supplied \code{tunpara} sequence. Typical usage is to have the program compute its own \code{tunpara} sequence based on \code{tunpara.num} and \code{tunpara.min.ratio}. Supplying a value of \code{tunpara} overrides this.
#' @param evatype the information criterion that will be used in the selection of optimal tuning parameter.
#' @param learnrate the learning rate in solving the optimization problem.
#' @param maxit specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{maxit = 5e3}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param threshold a threshold value. The deafult value is \code{1e-5}. If the absolute value of the each component of coefficients is smaller than \code{threshold}, it will be converted into zero.
#' @param trace a logical value. The default specification is \code{TRUE}, indicating that several useful information about the current process of fitting will be displayed.
#'
#' @return The fitted results are returned (a list). The estimates of regression coefficients is \code{res}.
#'
#' @examples
#' ### ==== An example for doing online variable selection via the promotion time cure model ==== ###
#'
#' ## ---- generate the simulated dataset ---- ##
#' sdata.full  <- sdata.PTC.High(N=500)
#' yobs  <- sdata.full$sdata$yobs
#' delta <- sdata.full$sdata$delta
#' X     <- as.matrix(sdata.full$sdata$X)
#' rm(sdata.full)
#'
#' ## ---- fit the promotion time cure model ---- ##
#' # estimator of the regression coefficients
#' fit.PTC <- PTC.NPMLE.ALasso.fit(yobs,delta,X)
#' print(fit.PTC$res)
#' # estimator of the baseline distribution at certain points
#' tm <- c(0.1,0.3,0.5,0.7,0.9)
#' Est.Lam0t <- PTC.NPMLE.Lam(tm=tm,bet=fit.PTC$res[,1],yobs=yobs,delta=delta,X=X)
#' print(data.frame(Est=Est.Lam0t,row.names=tm))
#'
#' @export PTC.NPMLE.ALasso.fit
PTC.NPMLE.ALasso.fit <- function(yobs,delta,X,LinkFunc=list(D0=exp,D1=exp),
                                 tunpara.num=50,tunpara.min.ratio=0.0001,tunpara=NULL,
                                 evatype=c("AIC","BIC")[2],
                                 learnrate=0.1,
                                 maxit=1e3,eps=1e-6, # the maximum number and tolerance of iteration in proximal gradient descent
                                 threshold=1e-5,     # the threshold for nonzero elements in lasso
                                 trace=TRUE){

  ### Preparations
  N <- length(yobs)
  XI <- cbind(1,X)
  pbet <- ncol(XI)

  ### calculate initial values for bet (no penalization ones)
  fit.nopen <- PTC.NPMLE.fit(yobs,delta,X,LinkFunc=LinkFunc,
                             maxit=maxit,eps=eps,simplify=TRUE)
  res.nopen <- fit.nopen$res
  bet.nopen <- res.nopen[,1]
  InfoM.nopen <- fit.nopen$InfoM
  weights  <- abs(1/bet.nopen[-1])

  ### obtain the basic alasso penalized parameter: using proximal gradient descent algorithm

  ## prepare candidate set of tuning parameters
  if(is.null(tunpara)){
    ScoreInfoM <- PTC.NPMLE.ScoreInfoM(bet=bet.nopen,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=TRUE)
    tunpara.max <-  max(abs((diag(ScoreInfoM$InfoM)*bet.nopen+ScoreInfoM$Score)[-1]*bet.nopen[-1]))
    tunpara <- exp(seq(log(tunpara.max*1.5), log(tunpara.max * tunpara.min.ratio),
                       length.out = tunpara.num))
    # tunpara <- sqrt(log((pbet))/N)*rev(tunscale) # to candidate set of tuning parameters
    # tunpara.num <- length(tunpara)
  }else{
    tunpara.num <- length(tunpara)
  }

  ## prepare boxes for puting results
  bet.path <- array(0,dim=c(2+pbet,tunpara.num),dimnames=list(c("IsConverge","EvaC",paste("bet",0:(pbet-1),sep="")),
                                                              paste("Tuning",1:tunpara.num,sep="")))

  ## iteratively solving SCAD for each tuning parameter
  method.solve <- "Taylor"
  bet.init <- c(bet.nopen[1],rep(0,pbet-1))
  if(trace==TRUE){cat("Fit the adaptive Lasso estimator: \n")}
  for(itunpara in 1:tunpara.num){ # itunpara <- 4
    if(trace==TRUE){cat(" - Current tuning parameter",itunpara,"/",length(tunpara),"(",tunpara[itunpara],") ... ")}

    # fit the lasso for current tuning parameter
    if(method.solve=="ProxGrad"){
      alasso.fit.c <- PTC.NPMLE.ALasso.ProxGrad(
        yobs=yobs,delta=delta,X=X,bet.init=bet.init,weights=weights,tunpara=tunpara[itunpara],LinkFunc=LinkFunc,
        learnrate=learnrate,maxit=maxit,eps=eps)
    }else if(method.solve=="Taylor"){
      alasso.fit.c <- PTC.NPMLE.ALasso.Taylor(
        yobs=yobs,delta=delta,X=X,bet.init=bet.init,weights=weights,tunpara=tunpara[itunpara],LinkFunc=LinkFunc,
        learnrate=learnrate,maxit=maxit,eps=eps)
    }
    bet.c <- alasso.fit.c$bet.alasso; round(bet.c,5)

    # calculate evaluating criterion
    EvaC.c <- as.numeric(-PTC.NPMLE.Beta.LogLik(bet.c,yobs,delta,XI,LinkFunc)+0.5*ifelse(evatype=="AIC",2,log(N))*sum(abs(bet.c[-1])>=threshold) )/N
    cat(EvaC.c,"\n")
    # combine and update the estimates
    bet.path[,itunpara] <- c(alasso.fit.c$convergence,EvaC.c,bet.c)
    bet.init <- bet.c
  }

  ## choose the best result according the chosen EvaC
  EvaC.best.idx <- which.min(bet.path['EvaC',])
  bet <- bet.path[-c(1,2),EvaC.best.idx]; tunpara.opt <- tunpara[EvaC.best.idx]
  EvaC <- bet.path[2,EvaC.best.idx]; convergence <- bet.path[1,EvaC.best.idx]==1
  round(bet,5)

  ### Obtain estimates of standard error
  ScoreInfoM <- PTC.NPMLE.ScoreInfoM(bet=bet,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=TRUE)
  Score <- ScoreInfoM$Score
  InfoM <- ScoreInfoM$InfoM
  # method 1:
  InfoM11 <- InfoM[bet!=0,bet!=0]
  InfoM12 <- InfoM[bet!=0,bet==0]
  InfoM22 <- InfoM[bet==0,bet==0]
  InfoM11.inv <- MASS::ginv(InfoM11)
  InfoM11.tilde.inv <- MASS::ginv(InfoM11 + tunpara.opt*diag(1/bet[bet!=0]^2))
  E.inv <- MASS::ginv(InfoM22 - t(InfoM12)%*%InfoM11.inv%*%InfoM12)
  covbet1 <- InfoM11.inv + (InfoM11.inv-InfoM11.tilde.inv)%*%InfoM12%*%E.inv%*%t(InfoM12)%*%(InfoM11.inv-InfoM11.tilde.inv)
  bet.se <- rep(NA,pbet); bet.se[bet!=0] <- sqrt(diag(covbet1)/N); bet.se
  # # method 2:
  # bet.se <- rep(NA,pbet); bet.se[bet!=0] <- sqrt(diag(MASS::ginv(InfoM[bet!=0,bet!=0]))/N); bet.se

  ### tidy the results: inference
  zvalue.bet <- bet/bet.se
  pvalue.bet <- 2*(1-pnorm(abs(zvalue.bet)))
  res <- data.frame(Est=bet,SE=bet.se,zvalue=zvalue.bet,pvalue=pvalue.bet,
                    row.names=c("Intercept",colnames(X)))

  ### output
  out <- list(
    convergence = convergence,
    bet.path = bet.path,
    res = res,
    Score = Score,
    InfoM = InfoM,
    res.nopen = res.nopen,
    bet.nopen = bet.nopen,
    InfoM.nopen = InfoM.nopen
  )
  return(out)

}


#==== Fit the adaptive Lasso (one time, for a pre-specified tuning parameter, using Taylor expansion) ====#
#' @title Variable selection in promotion time cure model (via adaptive Lasso, only one time for a specific tuning parameter, using Taylor expansion)
#'
#' @description Conduct variable selection using promotion time cure model based on the nonparametric maximum likelihood estimation (NPMLE) procedure and adaptive Lasso penalty (using Taylor expansion).
#' Note that it is a auxiliary function and only solve the optimization problem for a specific tuning parameter.
#'
#' @aliases PTC.NPMLE.ALasso.Taylor
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet.init the initial value for the regression coefficients.
#' @param weights the weight vector for the adaptive Lasso.
#' @param tunpara the current tuning parameter (scaler).
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' @param learnrate the learning rate in solving the optimization problem.
#' @param maxit specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{maxit = 5e3}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#'
#' @export PTC.NPMLE.ALasso.Taylor
PTC.NPMLE.ALasso.Taylor <- function(yobs,delta,X,bet.init,weights=weights,tunpara,
                                    LinkFunc=list(D0=exp,D1=exp),
                                    learnrate=0.1,maxit=5e3,eps=1e-6){

  ### Preparations
  pbet <- length(bet.init)
  tunpara.weights <- c(0,tunpara*weights)

  ### fit the model based on the current tuning parameter
  numit <- 1; bet.old <- bet.init
  repeat{

    ScoreInfoM <- PTC.NPMLE.ScoreInfoM(bet=bet.old,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=TRUE)
    XXXtYYY <- ScoreInfoM$InfoM%*%bet.old+ScoreInfoM$Score
    # XXX <- chol(ScoreInfoM$InfoM)
    # YYY <- as.vector(MASS::ginv(t(XXX)) %*% (ScoreInfoM$InfoM%*%bet.old+ScoreInfoM$Score))

    numit.inner <- 1
    bet.inner.old <- bet.old
    repeat{

      dev <- learnrate * ( ScoreInfoM$InfoM%*%bet.inner.old - XXXtYYY )
      bet.inner.temp <- bet.inner.old - as.vector(dev)
      bet.inner <- Soft.Threshold(bet.inner.temp,learnrate*tunpara.weights)

      # judge the convergence (inner)
      if( max(abs(bet.inner-bet.inner.old))>eps & numit.inner<maxit ){
        bet.inner.old <- bet.inner
        numit.inner <- numit.inner + 1
      }else{
        bet <- bet.inner
        break
      }

    }

    # judge the convergence (inner)
    if( max(abs(bet-bet.old))>eps & numit<maxit ){
      bet.old <- bet
      numit <- numit + 1
    }else{
      break
    }
  }
  convergence <- (numit<maxit)

  ### output
  out <- list(
    convergence = convergence,
    bet.alasso = bet
  )
  return(out)

}


#==== Fit the adaptive Lasso (one time, for a pre-specified tuning parameter, using proximal gradient descent) ====#
#' @title Variable selection in promotion time cure model (via adaptive Lasso, only one time for a specific tuning parameter, using proximal gradient descent)
#'
#' @description Conduct variable selection using promotion time cure model based on the nonparametric maximum likelihood estimation (NPMLE) procedure and adaptive Lasso penalty (using proximal gradient descent).
#' Note that it is a auxiliary function and only solve the optimization problem for a specific tuning parameter.
#'
#' @aliases PTC.NPMLE.ALasso.ProxGrad
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet.init the initial value for the regression coefficients.
#' @param weights the weight vector for the adaptive Lasso.
#' @param tunpara the current tuning parameter (scaler).
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' @param learnrate the learning rate in solving the optimization problem.
#' @param maxit specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{maxit = 5e3}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#'
#' @export PTC.NPMLE.ALasso.ProxGrad
PTC.NPMLE.ALasso.ProxGrad <- function(yobs,delta,X,bet.init,weights,tunpara,
                                      LinkFunc=list(D0=exp,D1=exp),
                                      learnrate=0.1,maxit=5e3,eps=1e-6){

  ### Preparations

  ### fit the model based on the current tuning parameter
  numit <- 1;
  bet.old <- bet.init
  repeat{
    Score <- PTC.NPMLE.ScoreInfoM(bet=bet.old,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=FALSE)$Score
    dev <- learnrate * Score
    bet.temp <- bet.old + as.vector(dev)
    bet.alasso <- c(bet.temp[1],Soft.Threshold(theta=bet.temp[-1],lam=weights*tunpara*learnrate))
    # print(round(bet.alasso[1:8],6))
    # print(max(abs(bet.alasso-bet.old)))
    if( max(abs(bet.alasso-bet.old))>eps & numit<maxit ){
      bet.old <- bet.alasso
      numit <- numit + 1
    }else{
      break
    }
  }
  convergence <- (numit<maxit)

  ### output
  out <- list(
    convergence = convergence,
    bet.alasso = bet.alasso
  )
  return(out)

}


#==== The soft-thresholding function ====#
#' @title The soft-thresholding function
#'
#' @description The soft-thresholding function.
#'
#' @aliases Soft.Threshold
#'
#' @param theta the parameter that we will calculate at.
#' @param lam the inputted tuning parameter
#'
#' @export Soft.Threshold
Soft.Threshold <- function(theta,lam){
  # the solution to LASSO penalty under the simplest situation: 2^{-1}(z-theta)^2+penalty
  # this is the equation (2.6) in Fan and Li (2001)
  res <- sign(theta)*pmax(abs(theta)-lam,0)
  return(res)
}


