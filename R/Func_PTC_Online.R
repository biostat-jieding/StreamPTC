



#==========================================================================#
# Fit the promotion time cure model with massive dataset
# -- based on: NPMLE from Portier et al. 2017
# -- a Renewable approach (basic estimation or with variable selection procedure)
#==========================================================================#


#==== Function for fitting PTC with NPMLE (one batch, Renewable approach) ====#
#' @title Online estimation in promotion time cure model
#' based on the classical nonparametric maximum likelihood estimator
#'
#' @description Fit online promotion time cure model based on the nonparametric maximum likelihood estimation (NPMLE) procedure.
#' In other words, this function can analyze streaming survival data with a cure fraction.
#' More details about the classical NPMLE method can be found in Portier et al. (2017).
#'
#' @aliases PTC.NPMLE.Renew.Batch
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param initial a logical value. The default specification is \code{TRUE}, indicating that the current data batch is the first batch and no historical summary statistics is available.
#' @param tm The time points that the nonparametric baseline distribution will be estimated at. The default specification is \code{TRUE}.
#' @param prev indicates the historical summary data. The default specification is \code{NULL} (with \code{initial=TRUE}). Otherwise, it is a list with the following five elements:
#'   \code{bet} is the historical estimated result of regression coefficients;
#'   \code{InfoM} is the historical estimated result of the information matrix;
#'   \code{N} is the sample of historical raw data;
#'   \code{Lamt} is the historical estimated result of the nonparametric baseline distribution at tm (if \code{tm=NULL}, we can set \code{Lamt=NULL} directly);
#'   \code{LamtInfoM} is the historical estimated result of the information matrix for the nonparametric baseline distribution at tm (if \code{tm=NULL}, we can set \code{LamtInfoM=NULL} directly.)
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#' @param maxit specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{maxit = 5e3}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param simplify a logical value. The default specification is \code{TRUE}, indicating that the inversion at each Newton iteration step will be avoided (using the inversion of the information matrix from historical data directly).
#'
#' @return The fitted results are returned (a list). The estimates of regression coefficients is \code{res}.
#'
#' @references PORTIER, F., EL GHOUCH, A. and VAN KEILEGOM, I. (2017). Efficiency and bootstrap in the promotion time cure model. Bernoulli 23 3437–3468.
#'
#' @examples
#' ### ==== An example for fitting the online promotion time cure model ==== ###
#'
#' ## ---- generate the dataset (full) ---- ##
#' N <- 10000
#' sdata.full  <- sdata.PTC(N=N)
#' yobs  <- sdata.full$sdata$yobs
#' delta <- sdata.full$sdata$delta
#' X     <- as.matrix(sdata.full$sdata$X)
#' rm(sdata.full)
#'
#' ## ---- fit the promotion time cure model in an online manner ---- ##
#' # prepare basic elements
#' tm <- c(0.1,0.3,0.5,0.7,0.9)
#' B <- 40
#' Nb <- N/B
#' batch <- ceiling((1:N) / Nb)
#' Res.Online <- list(
#'   res = array(0,dim=c(B,3,2),
#'               dimnames=list(paste("batch",1:B,sep=""),paste("bet",0:2,sep=""),c("EST","SE"))),
#'   res.Lamt = array(0,dim=c(B,length(tm)),
#'                    dimnames=list(paste("batch",1:B,sep=""),paste("tm",1:length(tm),sep="")))
#' )
#' # the online procedure
#' for(b in 1:B){  # b <- 1
#'   cat("Batch",b,"...\n")
#'
#'   # preparation: the batch idx / the previous elements
#'   idxb <- batch==b
#'   if(b==1){ prevb <- NULL }else{
#'     prevb <- list(
#'       bet       = fitb$res[,1],
#'       InfoM     = fitb$InfoM,
#'       Lamt      = fitb$Lamt,
#'       LamtInfoM = fitb$LamtInfoM,
#'       N         = fitb$N
#'     )
#'   }
#'
#'   # fit the current data batch (with current data the historical statistics)
#'   fitb <- PTC.NPMLE.Renew.Batch(
#'     yobs=yobs[idxb],delta=delta[idxb],X=as.matrix(X[idxb,,drop=FALSE]),
#'     initial=(b==1),tm=tm,prev=prevb)
#'   Res.Online$res[b,,1] <- fitb$res[,1]
#'   Res.Online$res[b,,2] <- fitb$res[,2]
#'   Res.Online$res.Lamt[b,] <- fitb$res.Lamt[,1]
#'
#' }
#' # present the fitted results
#' print(Res.Online$res)
#' print(Res.Online$res.Lamt)
#'
#' @export PTC.NPMLE.Renew.Batch
PTC.NPMLE.Renew.Batch <- function(yobs,delta,X,initial=TRUE,tm=NULL,prev=NULL,
                                  LinkFunc=list(D0=exp,D1=exp),maxit=1e3,eps=1e-6,simplify=TRUE){

  ### prepare
  XI <- cbind(1,X)
  pbet <- ncol(XI)
  N <- length(yobs)

  ### calculate the updated estimator
  if(initial==TRUE){

    ### fit the model using classical method
    fit.init <- PTC.NPMLE.fit(yobs,delta,X,LinkFunc=LinkFunc,
                              maxit=maxit,eps=eps,simplify=simplify)
    # for bet
    bet <- fit.init$res[,1]
    bet.se <- fit.init$res[,2]
    Score <- fit.init$Score
    InfoM <- fit.init$InfoM
    # for Lam
    if(!is.null(tm)){
      Lamt <- PTC.NPMLE.Lam(tm=tm,bet=bet,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc)
      LamtInfoM <- NA
      Lamt.se <- NA
    }

  }else{

    ### do repeation (for better InfoM)
    bet.prev       <- prev$bet
    InfoM.prev     <- prev$InfoM
    N.prev         <- prev$N
    if(!is.null(tm)){
      Lamt.prev      <- prev$Lamt
      LamtInfoM.prev <- prev$LamtInfoM
    }

    # updated estimator based on previous summary information and current data batch
    if(simplify==TRUE){

      # based on newton method: similar to Luo and Song (2020)
      solve.InfoM <- solve(
        N*PTC.NPMLE.ScoreInfoM(bet=bet.prev,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=FALSE,IsInfoM=TRUE)$InfoM +
          N.prev*InfoM.prev)
      numit <- 1
      bet.old <- bet.prev
      repeat{
        # ---- 1 ---- #
        Score <- PTC.NPMLE.ScoreInfoM(bet=bet.old,yob=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=FALSE)$Score
        dev <- solve.InfoM%*%(N.prev*InfoM.prev%*%(bet.prev-bet.old)+N*Score)
        # ---- 2 ---- #
        # InfoMScore <- PTC.NPMLE.ScoreInfoM(bet=bet.old,yob=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=TRUE)
        # Score <- InfoMScore$Score
        # InfoM <- InfoMScore$InfoM
        # dev <- solve(N*InfoM+N.prev*InfoM.prev)%*%(N.prev*InfoM.prev%*%(bet.prev-bet.old)+N*Score)
        # ---- 3 ---- #
        bet <- bet.old + as.vector(dev)
        if( max(abs(bet-bet.old))>eps & numit<maxit ){
          bet.old <- bet
          numit <- numit + 1
        }else{
          break
        }
      }

    }else{

      # based on the augmented likelihood directly
      sol <- stats::optim(par=rep(0,length(bet.prev)),fn=PTC.NPMLE.Renew.Beta.LogLik,
                          control=list(maxit=maxit,fnscale=-1,reltol=eps),
                          yobs=yobs,delta=delta,XI=XI,
                          bet.prev=bet.prev,InfoM.prev=InfoM.prev,ScoreAppro.prev=0,N.prev=N.prev,LinkFunc=LinkFunc)
      bet <- sol$par

    }

    ### calculate SEs for bet using explicit formula !!!
    # current elements
    InfoMScore.curr <- PTC.NPMLE.ScoreInfoM(bet=bet,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=TRUE)
    InfoM.curr      <- InfoMScore.curr$InfoM
    # update various elements
    InfoM      <- (N.prev/(N.prev+N))*InfoM.prev      + (N/(N.prev+N))*InfoM.curr
    # calculate SEs for bet
    bet.se <- sqrt(diag(solve(InfoM))/(N.prev+N))

    ### calculate SEs for Lamt using explicit formula !!!
    if(!is.null(tm)){
      # current elements
      Lamt.curr       <- PTC.NPMLE.Lam(tm=tm,bet=bet,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc)
      LamtInfoM.curr  <- NA
      # update various elements
      Lamt       <- (N.prev/(N.prev+N))*Lamt.prev       + (N/(N.prev+N))*Lamt.curr
      LamtInfoM  <- (N.prev/(N.prev+N))*LamtInfoM.prev  + (N/(N.prev+N))*LamtInfoM.curr
      # calculate SEs for bet and Lamt
      Lamt.se <- sqrt((1/LamtInfoM)/(N.prev+N))
    }

  }

  # summary information for this current step
  # for bet
  zvalue.bet <- bet/bet.se
  pvalue.bet <- 2*(1-pnorm(abs(zvalue.bet)))
  res <- data.frame(Est=bet,SE=bet.se,zvalue=zvalue.bet,pvalue=pvalue.bet,
                    row.names=c("(Intercept)",colnames(X)))
  # for Lam
  if(!is.null(tm)){
    zvalue.Lamt <- Lamt/Lamt.se
    pvalue.Lamt <- 2*(1-pnorm(abs(zvalue.Lamt)))
    res.Lamt    <- data.frame(Est=Lamt,SE=Lamt.se,zvalue=zvalue.Lamt,pvalue=pvalue.Lamt,
                              row.names=paste("tm=",tm,sep=""))
  }else{
    res.Lamt <- Lamt <- LamtInfoM <- NULL
  }


  ### output: coef matrix + next step's summary information ###
  out <- list(
    res = res,
    InfoM = InfoM,
    res.Lamt = res.Lamt,
    Lamt = Lamt,
    LamtInfoM = LamtInfoM,
    N = ifelse(initial==TRUE,N,N+N.prev)
  )
  return(out)


}


#==== High-dimensional Variable Selection in PTC with NPMLE (one batch, Renewable, ALasso) ====#
#' @title Online variable selection in promotion time cure model
#' based on the classical nonparametric maximum likelihood estimator and adaptive Lasso penalty
#'
#' @description Conduct online variable selection using promotion time cure model based on the nonparametric maximum likelihood estimation (NPMLE) procedure and adaptive Lasso penalty.
#' In other words, this function can analyze streaming survival data with a cure fraction and selection the underlying useful variables.
#'
#' @aliases PTC.NPMLE.Renew.ALasso.Batch
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param initial a logical value. The default specification is \code{TRUE}, indicating that the current data batch is the first batch and no historical summary statistics is available.
#' @param tm The time points that the nonparametric baseline distribution will be estimated at. The default specification is \code{TRUE}.
#' @param prev indicates the historical summary data. The default specification is \code{NULL} (with \code{initial=TRUE}). Otherwise, it is a list with the following five elements:
#'   \code{bet} is the historical estimated result of regression coefficients;
#'   \code{InfoM} is the historical estimated result of the information matrix;
#'   \code{N} is the sample of historical raw data;
#'   \code{Lamt} is the historical estimated result of the nonparametric baseline distribution at tm (if \code{tm=NULL}, we can set \code{Lamt=NULL} directly);
#'   \code{LamtInfoM} is the historical estimated result of the information matrix for the nonparametric baseline distribution at tm (if \code{tm=NULL}, we can set \code{LamtInfoM=NULL} directly.)
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
#' @param Include a logical value. The default specification is \code{TRUE}, indicating that the recovering term will be included.
#'
#' @return The fitted results are returned (a list). The estimates of regression coefficients is \code{res}.
#'
#' @examples
#' ### ==== An example for doing online variable selection via the promotion time cure model ==== ###
#'
#' ## ---- generate the dataset (full) ---- ##
#' N <- 4000
#' pX <- 10
#' sdata.full  <- sdata.PTC.High(N=N,pX=pX)
#' yobs  <- sdata.full$sdata$yobs
#' delta <- sdata.full$sdata$delta
#' X     <- as.matrix(sdata.full$sdata$X)
#' rm(sdata.full)
#'
#' ## ---- fit the promotion time cure model in an online manner ---- ##
#' # prepare basic elements
#' tm <- c(0.1,0.3,0.5,0.7,0.9)
#' B <- 20
#' Nb <- N/B
#' batch <- ceiling((1:N) / Nb)
#' Res.Online_ALasso <- list(
#'   res = array(0,dim=c(B,pX+1,2),
#'               dimnames=list(paste("batch",1:B,sep=""),paste("bet",0:pX,sep=""),c("EST","SE"))),
#'   res.Lamt = array(0,dim=c(B,length(tm)),
#'                    dimnames=list(paste("batch",1:B,sep=""),paste("tm",1:length(tm),sep="")))
#' )
#' # the online procedure
#' for(b in 1:B){  # b <- 1
#'   cat("Batch",b,"...\n")
#'
#'   # preparation: the batch idx / the previous elements
#'   idxb <- batch==b
#'   if(b==1){ prevb <- NULL }else{
#'     prevb <- list(
#'       bet         = fitb$res[,1],
#'       Score       = fitb$Score,
#'       InfoM       = fitb$InfoM,
#'       ScoreAppro  = fitb$ScoreAppro,
#'       bet.nopen   = fitb$bet.nopen,
#'       InfoM.nopen = fitb$InfoM.nopen,
#'       Lamt        = fitb$Lamt,
#'       LamtInfoM   = fitb$LamtInfoM,
#'       N           = fitb$N
#'     )
#'   }
#'
#'   # fit the current data batch (with current data the historical statistics)
#'   fitb <- PTC.NPMLE.Renew.ALasso.Batch(
#'     yobs=yobs[idxb],delta=delta[idxb],X=as.matrix(X[idxb,,drop=FALSE]),
#'     initial=(b==1),tm=tm,prev=prevb,Include=TRUE)
#'   Res.Online_ALasso$res[b,,1] <- fitb$res[,1]
#'   Res.Online_ALasso$res[b,,2] <- fitb$res[,2]
#'   Res.Online_ALasso$res.Lamt[b,] <- fitb$res.Lamt[,1]
#'
#' }
#' # present the fitted results
#' print(Res.Online_ALasso$res)
#' print(Res.Online_ALasso$res.Lamt)
#'
#' @export PTC.NPMLE.Renew.ALasso.Batch
PTC.NPMLE.Renew.ALasso.Batch <- function(yobs,delta,X,
                                         initial=TRUE,tm=NULL,prev=NULL,
                                         LinkFunc=list(D0=exp,D1=exp),
                                         tunpara.num=50,
                                         tunpara.min.ratio=0.0001,
                                         tunpara=NULL,
                                         evatype=c("AIC","BIC")[2],
                                         learnrate=c(0.1),
                                         maxit=1e3,
                                         eps=1e-6,
                                         threshold=1e-5,
                                         trace=TRUE,
                                         Include=TRUE # Include B0 to recover zeros
){

  ### calculate the updated estimator according to whether it is the first batch
  if(initial==TRUE){

    ### fit the model using classical method
    fit.init <- PTC.NPMLE.ALasso.fit(
      yobs=yobs,delta=delta,X=X,
      LinkFunc=LinkFunc,
      tunpara.num=tunpara.num,
      tunpara.min.ratio=tunpara.min.ratio,
      tunpara=tunpara,
      evatype=evatype,
      learnrate=learnrate,
      maxit=maxit,eps=eps,threshold=threshold,
      trace=trace
    )
    # for bet
    bet.alasso <- fit.init$res[,1]
    bet.alasso.se <- fit.init$res[,2]
    Score <- fit.init$Score
    InfoM <- fit.init$InfoM
    ScoreAppro <- Score
    # for nopen
    bet.nopen <- fit.init$bet.nopen
    InfoM.nopen <- fit.init$InfoM.nopen
    # for Lam
    Lamt <- PTC.NPMLE.Lam(tm=tm,bet=bet.alasso,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc)
    LamtInfoM <- NA
    Lamt.se <- NA
    # other
    N.prev <- 0
    N <- length(yobs)
    bet.path <- fit.init$bet.path

  }else if(initial==FALSE){

    ### do preparation for basic information
    XI <- cbind(1,X)
    pbet <- ncol(XI)
    N <- length(yobs)

    ### extract historical elements
    bet.alasso.prev  <- prev$bet
    Score.prev       <- prev$Score
    InfoM.prev       <- prev$InfoM
    ScoreAppro.prev  <- prev$ScoreAppro
    bet.nopen.prev   <- prev$bet.nopen
    InfoM.nopen.prev <- prev$InfoM.nopen
    Lamt.prev        <- prev$Lamt
    LamtInfoM.prev   <- prev$LamtInfoM
    N.prev           <- prev$N

    ### update the regression coefficients without considering penalty
    fit.nopen<- PTC.NPMLE.Renew.Batch(
      yobs=yobs,delta=delta,X=X,initial=FALSE,tm=0,
      prev=list(bet=bet.nopen.prev,InfoM=InfoM.nopen.prev,Lamt=0,LamtInfoM=0,N=N.prev),
      LinkFunc=LinkFunc,maxit=maxit,eps=eps)
    bet.nopen <- fit.nopen$res[,1]
    InfoM.nopen <- fit.nopen$InfoM
    weights <- abs(1/bet.nopen[-1])

    ### update basic adaptive lasso estimator: based on proximal gradient method

    ## prepare candidate set of tuning parameters and boxes for puting results
    if(is.null(tunpara)){
      ScoreInfoM <- PTC.NPMLE.ScoreInfoM(bet=bet.nopen,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=TRUE)
      InfoM.cc <- (N/(N.prev+N))*ScoreInfoM$InfoM + (N.prev/(N.prev+N))*InfoM.prev
      if(Include==TRUE){
        Score.cc <- (N.prev*ScoreAppro.prev+
                       N.prev*InfoM.prev%*%(bet.alasso.prev-bet.nopen)+
                       N*ScoreInfoM$Score)/(N+N.prev)
      }else{
        Score.cc <- (N.prev*InfoM.prev%*%(bet.alasso.prev-bet.nopen)+
                       N*ScoreInfoM$Score)/(N+N.prev)
      }
      tunpara.max <- max(abs((diag(InfoM.cc)*bet.nopen+Score.cc)[-1]*bet.nopen[-1]))
      tunpara <- exp(seq(log(tunpara.max*1.5), log(tunpara.max * tunpara.min.ratio),
                         length.out = tunpara.num))
      # tunpara <- sqrt(log((pbet))/(N+N.prev))*rev(tunscale)
      # tunpara.num <- length(tunpara)
    }else{
      tunpara.num <- length(tunpara)
    }


    ## prepare boxes for puting results
    bet.path <- array(0,dim=c(2+pbet,tunpara.num),dimnames=list(c("IsConverge","EvaC",paste("bet",0:(pbet-1),sep="")),
                                                                paste("Tuning",1:tunpara.num,sep="")))

    ## iteratively solving SCAD for each tuning parameter
    method.solve <- "Taylor"
    bet.init <- c(bet.alasso.prev[1],rep(0,pbet-1))
    if(trace==TRUE){cat("Fit the adaptive lasso estimator: \n")}
    for(itunpara in 1:tunpara.num){ # itunpara <- 3
      if(trace==TRUE){cat(" - Current tuning parameter",itunpara,"/",length(tunpara),"(",tunpara[itunpara],") ... ")}

      tunpara.weights <- c(0,tunpara[itunpara]*weights)

      # fit the model based on the current tuning parameter
      if(method.solve=="ProxGrad"){

        numit <- 1; bet.old <- bet.init
        repeat{
          Score <- PTC.NPMLE.ScoreInfoM(bet=bet.old,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=FALSE)$Score
          if(Include==TRUE){
            dev <- learnrate * (N.prev*ScoreAppro.prev+N.prev*InfoM.prev%*%(bet.alasso.prev-bet.old)+N*Score)/(N+N.prev)
          }else{
            dev <- learnrate * (N.prev*InfoM.prev%*%(bet.alasso.prev-bet.old)+N*Score)/(N+N.prev)
          }
          bet.temp <- bet.old + as.vector(dev)
          bet.c <- Soft.Threshold(theta=bet.temp,lam=tunpara.weights*learnrate)
          if( max(abs(bet.c-bet.old))>eps & numit<maxit ){
            bet.old <- bet.c
            numit <- numit + 1
          }else{
            break
          }
        }

      }else if(method.solve=="Taylor"){

        numit <- 1; bet.old <- bet.init
        repeat{

          ScoreInfoM <- PTC.NPMLE.ScoreInfoM(bet=bet.old,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=TRUE)
          InfoM.cc <- (N/(N.prev+N))*ScoreInfoM$InfoM + (N.prev/(N.prev+N))*InfoM.prev
          if(Include==TRUE){
            Score.cc <- (N.prev*ScoreAppro.prev+
                           N.prev*InfoM.prev%*%(bet.alasso.prev-bet.old)+
                           N*ScoreInfoM$Score)/(N+N.prev)
          }else{
            Score.cc <- (N.prev*InfoM.prev%*%(bet.alasso.prev-bet.old)+
                           N*ScoreInfoM$Score)/(N+N.prev)
          }
          XXXtYYY <- InfoM.cc%*%bet.old+Score.cc

          numit.inner <- 1
          bet.inner.old <- bet.old
          repeat{
            # use Proximal Gradient Descent Here
            dev <- learnrate * ( InfoM.cc%*%bet.inner.old - XXXtYYY )
            bet.inner.temp <- bet.inner.old - as.vector(dev)
            bet.inner <- Soft.Threshold(bet.inner.temp,learnrate*tunpara.weights)
            # judge the convergence (inner)
            if( max(abs(bet.inner-bet.inner.old))>eps & numit.inner<maxit ){
              bet.inner.old <- bet.inner
              numit.inner <- numit.inner + 1
            }else{
              bet.c <- bet.inner
              break
            }
          }

          # judge the convergence (inner)
          if( max(abs(bet.c-bet.old))>eps & numit<maxit ){
            bet.old <- bet.c
            numit <- numit + 1
          }else{
            break
          }
        }

      }
      # calculate evaluating criterion
      EvaC.c <- as.numeric(
        -PTC.NPMLE.Renew.Beta.LogLik(bet.c,yobs,delta,XI,bet.alasso.prev,InfoM.prev,ScoreAppro.prev,N.prev,LinkFunc)+
          0.5*ifelse(evatype=="AIC",2,log(N+N.prev))*sum(abs(bet.c[-1])>=threshold) )/(N+N.prev)
      cat(EvaC.c,"\n")
      # combine and update the estimates
      bet.path[,itunpara] <- c(numit<maxit,EvaC.c,bet.c)
      bet.init <- bet.c
    }

    ## choose the best result according the chosen EvaC
    EvaC.best.idx <- which.min(bet.path['EvaC',])
    bet.alasso <- bet.path[-c(1,2),EvaC.best.idx]; tunpara.opt <- tunpara[EvaC.best.idx]
    EvaC <- bet.path[2,EvaC.best.idx]; convergence <- bet.path[1,EvaC.best.idx]==1
    round(bet.alasso,5)

    ### calculating desparsified lasso estimator
    # extract basic elements
    InfoMScore.curr <- PTC.NPMLE.ScoreInfoM(bet=bet.alasso,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc,IsScore=TRUE,IsInfoM=TRUE)
    Score.curr    <- InfoMScore.curr$Score
    InfoM.curr    <- InfoMScore.curr$InfoM
    # update information matrix and score to renewable ones
    Score      <- (N.prev/(N.prev+N))*Score.prev    + (N/(N.prev+N))*Score.curr
    InfoM      <- (N.prev/(N.prev+N))*InfoM.prev    + (N/(N.prev+N))*InfoM.curr
    ScoreAppro <- as.vector( (N.prev/(N.prev+N))*(ScoreAppro.prev+InfoM.prev%*%(bet.alasso.prev-bet.alasso)) + (N/(N.prev+N))*Score.curr )
    # Obtain estimates of standard error
    # --- method 1:
    InfoM11 <- InfoM[bet.alasso!=0,bet.alasso!=0]
    InfoM12 <- InfoM[bet.alasso!=0,bet.alasso==0]
    InfoM22 <- InfoM[bet.alasso==0,bet.alasso==0]
    InfoM11.inv <- MASS::ginv(InfoM11)
    if(any(bet.alasso==0)){
      InfoM11.tilde.inv <- MASS::ginv(InfoM11 + tunpara.opt*diag(1/bet.alasso[bet.alasso!=0]^2))
      E.inv <- MASS::ginv(InfoM22 - t(InfoM12)%*%InfoM11.inv%*%InfoM12)
      covbet1 <- InfoM11.inv + (InfoM11.inv-InfoM11.tilde.inv)%*%InfoM12%*%E.inv%*%t(InfoM12)%*%(InfoM11.inv-InfoM11.tilde.inv)
    }else{
      covbet1 <- InfoM11.inv
    }
    bet.alasso.se <- rep(NA,pbet); bet.alasso.se[bet.alasso!=0] <- sqrt(diag(covbet1)/(N+N.prev)); bet.alasso.se
    # # --- method 2:
    # bet.alasso.se <- rep(NA,pbet); bet.alasso.se[bet.alasso!=0] <- sqrt(diag(MASS::ginv(InfoM[bet.alasso!=0,bet.alasso!=0]))/(N+N.prev))

    ### calculate SEs for Lamt using explicit formula !!!
    # current elements
    Lamt.curr       <- PTC.NPMLE.Lam(tm=tm,bet=bet.alasso,yobs=yobs,delta=delta,X=X,LinkFunc=LinkFunc)
    LamtInfoM.curr  <- NA
    # update various elements
    Lamt       <- (N.prev/(N.prev+N))*Lamt.prev       + (N/(N.prev+N))*Lamt.curr
    LamtInfoM  <- (N.prev/(N.prev+N))*LamtInfoM.prev  + (N/(N.prev+N))*LamtInfoM.curr
    # calculate SEs for bet and Lamt
    Lamt.se <- sqrt((1/LamtInfoM)/(N.prev+N))

  }

  # summary information for this current step
  # for bet
  zvalue.bet.alasso <- bet.alasso/bet.alasso.se
  pvalue.bet.alasso <- 2*(1-pnorm(abs(zvalue.bet.alasso)))
  res <- data.frame(Est=bet.alasso,SE=bet.alasso.se,zvalue=zvalue.bet.alasso,pvalue=pvalue.bet.alasso,
                    row.names=c("Intercept",colnames(X)))
  # for Lam
  zvalue.Lamt <- Lamt/Lamt.se
  pvalue.Lamt <- 2*(1-pnorm(abs(zvalue.Lamt)))
  res.Lamt    <- data.frame(Est=Lamt,SE=Lamt.se,zvalue=zvalue.Lamt,pvalue=pvalue.Lamt,
                            row.names=paste("tm=",tm,sep=""))

  ### output: coef matrix + next step's summary information ###
  out <- list(
    res = res,
    res.Lamt = res.Lamt,
    # for bet's Score and InfoM
    Score     = Score,
    InfoM     = InfoM,
    ScoreAppro=ScoreAppro,
    # for nopen bet
    bet.nopen = bet.nopen,
    InfoM.nopen = InfoM.nopen,
    # for lambda
    Lamt = Lamt,
    LamtInfoM = LamtInfoM,
    # other
    N = ifelse(initial==TRUE,N,N+N.prev),
    bet.path = bet.path
  )
  return(out)


}



#==== Online Loss function for bet (Profiled online log-likelihood function) ====#
#' @title Loss function for bet in NPMLE (Online version)
#'
#' @description  Online loss function for the promotion cure model (Profiled online log-likelihood function in NPMLE).
#' Note that it highly depends on the function \code{\link{PTC.NPMLE.Beta.LogLik}} in offline scenario.
#'
#' @aliases PTC.NPMLE.Renew.Beta.LogLik
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param XI a matrix of covariates with 1 as the first column, that is, \code{XI=cbind(1,X)}.
#' @param bet.prev the historical summary data of regression coefficients.
#' @param InfoM.prev the historical summary data of the information matrix (inverse variance-covariance matrix).
#' @param ScoreAppro.prev the historical summary data of a term derived from the online variable selection procedure (the recovering term).
#' @param N.prev the sample of historical raw data
#' @param LinkFunc the link function specified in the promotion time cure model. The default is the exponential function (The proportional hazards cure model).
#'
#' @export PTC.NPMLE.Renew.Beta.LogLik
PTC.NPMLE.Renew.Beta.LogLik <- function(bet,yobs,delta,XI,bet.prev,InfoM.prev,ScoreAppro.prev,N.prev,
                                        LinkFunc=list(D0=exp,D1=exp)){

  # calculate the loss function for beta (log form)
  val.log <- PTC.NPMLE.Beta.LogLik(bet=bet,yobs=yobs,delta=delta,XI=XI,LinkFunc=LinkFunc)
  val.log.online <-
    val.log -
    as.vector(t(bet-bet.prev)%*%(InfoM.prev%*%(bet-bet.prev)-2*ScoreAppro.prev)) * N.prev / 2


  # output
  return(val.log.online)

}





