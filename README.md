# StreamPTC
An R package for fitting promotion time cure (PTC) model in the streaming data environment.

## Four main R functions are included (with different scenarios for application):
- Fit PTC model in an offline manner under low-dimensional scenario
  - using *PTC.NPMLE.fit()*
- Fit PTC model in an offline manner under high-dimensional scenario (variable selection)
  - using *PTC.NPMLE.ALasso.fit()*
- Fit PTC model in an online manner under low-dimensional scenario
  - using *PTC.NPMLE.Renew.Batch()*
- Fit PTC model in an online manner under high-dimensional scenario (variable selection)
  - using *PTC.NPMLE.Renew.ALasso.Batch()*
  
## The corresponding examples using simulated datasets are shown below:

- *An example for fitting the promotion time cure model*
```R
## ---- generate the simulated dataset ---- ##
sdata.full  <- sdata.PTC(N=200)
yobs  <- sdata.full$sdata$yobs
delta <- sdata.full$sdata$delta
X     <- as.matrix(sdata.full$sdata$X)
rm(sdata.full)

## ---- fit the promotion time cure model ---- ##
# estimator of the regression coefficients
fit.PTC <- PTC.NPMLE.fit(yobs,delta,X)
print(fit.PTC$res)
# estimator of the baseline distribution at certain points
tm <- c(0.1,0.3,0.5,0.7,0.9)
Est.Lam0t <- PTC.NPMLE.Lam(tm=tm,bet=fit.PTC$res[,1],yobs=yobs,delta=delta,X=X)
print(data.frame(Est=Est.Lam0t,row.names=tm))
```

- *An example for doing online variable selection via the promotion time cure model*
```R
## ---- generate the simulated dataset ---- ##
sdata.full  <- sdata.PTC.High(N=500)
yobs  <- sdata.full$sdata$yobs
delta <- sdata.full$sdata$delta
X     <- as.matrix(sdata.full$sdata$X)
rm(sdata.full)

## ---- fit the promotion time cure model ---- ##
# estimator of the regression coefficients
fit.PTC <- PTC.NPMLE.ALasso.fit(yobs,delta,X)
print(fit.PTC$res)
# estimator of the baseline distribution at certain points
tm <- c(0.1,0.3,0.5,0.7,0.9)
Est.Lam0t <- PTC.NPMLE.Lam(tm=tm,bet=fit.PTC$res[,1],yobs=yobs,delta=delta,X=X)
print(data.frame(Est=Est.Lam0t,row.names=tm))
```

- *An example for fitting the online promotion time cure model*
```R
## ---- generate the dataset (full) ---- ##
N <- 10000
sdata.full  <- sdata.PTC(N=N)
yobs  <- sdata.full$sdata$yobs
delta <- sdata.full$sdata$delta
X     <- as.matrix(sdata.full$sdata$X)
rm(sdata.full)

## ---- fit the promotion time cure model in an online manner ---- ##
# prepare basic elements
tm <- c(0.1,0.3,0.5,0.7,0.9)
B <- 40
Nb <- N/B
batch <- ceiling((1:N) / Nb)
Res.Online <- list(
  res = array(0,dim=c(B,3,2),
              dimnames=list(paste("batch",1:B,sep=""),paste("bet",0:2,sep=""),c("EST","SE"))),
  res.Lamt = array(0,dim=c(B,length(tm)),
                   dimnames=list(paste("batch",1:B,sep=""),paste("tm",1:length(tm),sep="")))
)
# the online procedure
for(b in 1:B){  # b <- 1
  cat("Batch",b,"...\n")
  
  # preparation: the batch idx / the previous elements
  idxb <- batch==b
  if(b==1){ prevb <- NULL }else{
    prevb <- list(
      bet       = fitb$res[,1],
      InfoM     = fitb$InfoM,
      Lamt      = fitb$Lamt,
      LamtInfoM = fitb$LamtInfoM,
      N         = fitb$N
    )
  }
  
  # fit the current data batch (with current data the historical statistics)
  fitb <- PTC.NPMLE.Renew.Batch(
    yobs=yobs[idxb],delta=delta[idxb],X=as.matrix(X[idxb,,drop=FALSE]),
    initial=(b==1),tm=tm,prev=prevb)
  Res.Online$res[b,,1] <- fitb$res[,1]
  Res.Online$res[b,,2] <- fitb$res[,2]
  Res.Online$res.Lamt[b,] <- fitb$res.Lamt[,1]
  
}
# present the fitted results
print(Res.Online$res)
print(Res.Online$res.Lamt)
```

- *An example for doing online variable selection via the promotion time cure model*
```R
## ---- generate the dataset (full) ---- ##
N <- 4000
pX <- 10
sdata.full  <- sdata.PTC.High(N=N,pX=pX)
yobs  <- sdata.full$sdata$yobs
delta <- sdata.full$sdata$delta
X     <- as.matrix(sdata.full$sdata$X)
rm(sdata.full)

## ---- fit the promotion time cure model in an online manner ---- ##
# prepare basic elements
tm <- c(0.1,0.3,0.5,0.7,0.9)
B <- 20
Nb <- N/B
batch <- ceiling((1:N) / Nb)
Res.Online_ALasso <- list(
  res = array(0,dim=c(B,pX+1,2),
              dimnames=list(paste("batch",1:B,sep=""),paste("bet",0:pX,sep=""),c("EST","SE"))),
  res.Lamt = array(0,dim=c(B,length(tm)),
                   dimnames=list(paste("batch",1:B,sep=""),paste("tm",1:length(tm),sep="")))
)
# the online procedure
for(b in 1:B){  # b <- 1
  cat("Batch",b,"...\n")
  
  # preparation: the batch idx / the previous elements
  idxb <- batch==b
  if(b==1){ prevb <- NULL }else{
    prevb <- list(
      bet         = fitb$res[,1],
      Score       = fitb$Score,
      InfoM       = fitb$InfoM,
      ScoreAppro  = fitb$ScoreAppro,
      bet.nopen   = fitb$bet.nopen,
      InfoM.nopen = fitb$InfoM.nopen,
      Lamt        = fitb$Lamt,
      LamtInfoM   = fitb$LamtInfoM,
      N           = fitb$N
    )
  }
  
  # fit the current data batch (with current data the historical statistics)
  fitb <- PTC.NPMLE.Renew.ALasso.Batch(
    yobs=yobs[idxb],delta=delta[idxb],X=as.matrix(X[idxb,,drop=FALSE]),
    initial=(b==1),tm=tm,prev=prevb,Include=TRUE)
  Res.Online_ALasso$res[b,,1] <- fitb$res[,1] 
  Res.Online_ALasso$res[b,,2] <- fitb$res[,2]
  Res.Online_ALasso$res.Lamt[b,] <- fitb$res.Lamt[,1]
  
}
# present the fitted results
print(Res.Online_ALasso$res)
print(Res.Online_ALasso$res.Lamt)
```
