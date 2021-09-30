rm(list=ls())
setwd("~/Desktop/")

##1##    IMPORT & TRANSFORM DATA  ##
macro = ts(read.table("Macro M.csv", header = TRUE,sep = ";", 
                      skipNul = TRUE)[,-1], start = c(2001,01), frequency = 12)
google = ts(read.table("Google data.csv", header = TRUE,sep = ";", 
                       skipNul = TRUE)[,-1], start = c(2004,01), frequency = 12)
  #1a#   Seasonal Adjust    
install.packages("seasonal")
library(seasonal)

macro.sa = list()
for (i in (1:ncol(macro))) {
  if(i == 4 ){
    macro.sa[[i]] = macro[,i] #Government bond yield is already seasonally adjusted
    next
  } 
  if(i == 5 ){
    macro.sa[[i]] = macro[,i] #Unempoyment rate is already seasonally adjusted
    next
  }
  if(i == 6 ){
    macro.sa[[i]] = macro[,i] #Industrial Index is already seasonally adjusted
    next
  }
  macro.sa[[i]] = final(seas(macro[,i], na.action = na.x13, outlier = NULL))
}
macro.sa<- ts(matrix(unlist(macro.sa), ncol = 7, byrow = F), start = c(2001,01), frequency = 12)
  
  #1b# Log Difference
hstep = 1 
macrol = matrix(0, nrow= nrow(macro.sa)-hstep , ncol = ncol(macro.sa))
for (f in 1: ncol(macrol)){
  if(f == 4 ){
    macrol[,f] <- macro.sa[-(1:hstep),f] #Spread10Y-2Y 
    next
  }
  if(f == 5 ){
    macrol[,f] <- macro.sa[-(1:hstep),f] #Unemployment rate
    next
  }
  macrol[,f] <- (100/hstep)*diff(log(macro.sa[,f]), lag = hstep)
}

  #1c# Scaling to Zero mean, Unit variance
macrol <-ts(scale(macrol),start = c(2001,02), frequency = 12)
google <- scale(google)

  #1d# Matching two time series
macro.google <-  window(macrol,start = start(google))
google <-  window(macrol,start = start(google), end = end(macro.google))


##2## SETUP DEPENDENT AND EXPLANATORY VARIABLES   ##
dependent = vector("list", ncol(macro.google))
explanatoryX = vector("list", ncol(macro.google))   #Include lags of dependent var. and remain macro.var 
explanatoryZ = vector("list", ncol(macro.google))   #All Google.var corresponding to explanatory macro.var
google.prob = vector("list", ncol(google))          #All Google.Prob corresponding to explanatory macro.var
for (i in 1:ncol(macro.google)){
  dependent[[i]] = macro.google[3:nrow(macro.google),i]     #y_t
  explanatoryX[[i]] = cbind(
    macro.google[2:(nrow(macro.google)-1),i],               #y_t-1
    macro.google[1:(nrow(macro.google)-2),i],               #y_t-2
    macro.google[2:(nrow(macro.google)-1),-i]               #X_t-1
  )
  explanatoryZ[[i]]=google[3:nrow(macro.google),-i]         #Z_t 
  google.prob[[i]] = google[3:nrow(macro.google),-i]/100
}


##3## MODELS
  ##3A## With Google as Predictors
DM <- function(lambda,gamma,kappa,insample){
  dmG<- function(Y, x, z, insample, lambda , gamma, kappa){
    #3a# Linear regression
    linear <- lm(Y ~ x+z-1)
    m = ncol(x)
    variance.y<- summary(linear)$sigma
    
    #3b# Auto regression order 2 
    AR <- ar(Y, order.max = 2, method = "ols", demean = T)
    
    #3c# Time-Varying Parameters model 
    
    #3c.1Constructing the matrix to indicate which variables enter the model
    #For X
    models.whichX <- matrix(unlist(unique(combn(rep(c(0,1),m),m,simplify=FALSE))),(2^m),m,byrow=TRUE)
    models.whichX <- models.whichX[-(which(rowSums(models.whichX)==0)),]
    #For Z: Google variables conresponding to Macro. variables which enter the model 
    models.whichZ <- models.whichX[,-(1:2)]
    models.which <- cbind(models.whichX,models.whichZ)
    
    #3c.2# Initial matrices
    TT <- length(Y)
    K <- nrow (models.which)
    yhat <- matrix (rep(0,TT*K), ncol=K)
    theta<- vector ("list",K)
    sigma<- vector ("list",K)
    R <- vector ("list",K)
    theta0 <- rep(0,(ncol(x)+ncol(z)))
    
    eAll <- matrix (rep(0,TT*K), ncol=K)
    VAll <- rep(0,K)
    thetahat <-matrix (rep(0,K*(ncol(x)+ncol(z))),nrow=K)  
    Vtheta <- matrix (rep(0,K*(ncol(x)+ncol(z))),nrow=K)  
    yhat.Kalman <- matrix (rep(0,TT*K), ncol=K)
    predvar.Kalman <- rep(0,K)
    
    model.dim <- rep(0,K)
    for (k in 1:K) model.dim[k] <- sum(as.numeric(models.which[k,]))
    
    s2y <- variance.y # variance of Y
    s2x <- rep(0,(ncol(x)))
    s2z <- rep(0,(ncol(z)))
    for (a in 1:ncol(x)){
      for (b in 1:ncol(z)){
        s2x[a] <- var(x[,a])
        s2z[b] <- var(z[,b])
      }
    }
    s2X = c(s2x,s2z)
    for (k in 1:K){
      theta[[k]] <-linear$coefficients[as.logical( models.which[k,])]
      if (model.dim[k] == 1) 
        sigma[[k]] <- s2y/s2X[as.logical(models.which[k,])]
      else
        sigma[[k]] <- diag(s2y/s2X[as.logical(models.which[k,])])
    }
    #3c.3# TVP model
    rm.Kalman <-function (theta, sigma, V,lambda, X, Y, TT) {
      R = list()
      theta = list(); theta[[1]] <- theta0
      sigma = list(); sigma[[1]] <- sigma0
      update = list()
      e <- rep(0,TT)
      V <-  (rep(0,TT));V[1]<- V0
      A <-  (rep(0,TT))
      yhat <- rep(0,TT)
      for (tt in 2:(TT)){
        R[[tt]] <- sigma[[tt-1]]/lambda
        if(ncol(as.matrix(sigma0))==1){
          xtmp <- X[tt]
        } else{
          xtmp <- X[tt,]
        }
        update[[tt]] <- R[[tt]]%*% xtmp %*%  solve( V[tt] + t(xtmp) %*% R[[tt]] %*% xtmp)
        yhat[tt] <- t(xtmp) %*% theta[[tt-1]]
        predvar <- t(xtmp) %*% R[[tt]] %*% xtmp + V[tt]
        e[tt] <- Y[tt]-yhat[tt]
        theta[[tt]] <- theta[[tt-1]] +update[[tt]] %*% e[tt]
        sigma[[tt]] <- R[[tt]] - update[[tt]] %*% t(xtmp) %*% R[[tt]]
        V[tt+1] <- kappa*V[tt] + (1-kappa)*e[tt]^2
        ifelse(V[tt+1] <= 0,V[tt],V[tt+1])
      }
      list (theta=theta[[TT]], sigma=sigma[[TT]],e=e, yhat=yhat, predvar=predvar, V = V[TT],R=R[[TT]])
    }
    
    #3c.3# K-cases of TVP model
    Kalmanall = vector("list",K)
    X = cbind(x,z)
    for (k in 1:K){
      if (model.dim[k] == 1){
        sigma0 <- s2y/s2X[as.logical( models.which[k,])]
      } else {
        sigma0 <- diag(s2y/s2X[as.logical( models.which[k,])]) 
      }
      theta0 <- linear$coefficients[as.logical( models.which[k,])]
      V0 <- variance.y<- summary(linear)$sigma
      Xtmp <-  X[,as.logical( models.which[k,])]
      Kalmanall[[k]] <-rm.Kalman (theta[[k]], sigma[[k]], V,lambda, Xtmp, Y, TT)
      ###
      thetahat[k,as.logical( models.which[k,])] <- Kalmanall[[k]]$theta
      Vtheta[k,as.logical( models.which[k,])] <- diag(Kalmanall[[k]]$sigma)
      yhat.Kalman[,k] <- Kalmanall[[k]]$yhat
      predvar.Kalman[k] <- Kalmanall[[k]]$predvar
      eAll[,k] <- Kalmanall[[k]]$e
      VAll[k] <- Kalmanall[[k]]$V
    } 
    
    #3d# Model Probabilities
    eps=.001/nrow(models.which)
    model.update <- function (piold, gamma, eps,  y, yhat, predvar) {
      pipred <- (piold^gamma)+eps / sum((piold^gamma)+eps)
      logpyt <- -0.5*log(predvar) - 0.5*(y-yhat)^2/predvar # Update pi values
      pyt <- exp (logpyt)
      pinew <- pipred * pyt
      pinew <- pinew/sum(pinew)
      pinew <- pinew + eps
      pinew <- pinew/sum(pinew)
      
      list (pinew=as.vector(pinew), logpyt = logpyt)
    }
    
    pimat <- list()
    piOld <- rep (1/K, K)
    piMat <- matrix (rep(0,(TT)*K), ncol=K)
    logPL <- matrix (rep(0,(TT)*K), ncol=K)
    #3d.1# Update model probabilities
    for (tt in 1:(TT)){
      pimat[[tt]] <- model.update(piOld, gamma, eps, Y[tt], yhat.Kalman[tt,],predvar.Kalman)
      piMat[tt,] <- pimat[[tt]]$pinew
      logPL[tt,] <- pimat[[tt]]$logpyt
    }
    
    #3e# DMA modelling
    thetahat.ma <- matrix(0, nrow = TT, ncol = ncol(X))
    yhat.ma <- rep(0,TT)
    for (tt in 1:(TT)){
      thetahat.ma[tt,] <- t(piMat[tt,]) %*% thetahat
      yhat.ma[tt] <- t(piMat[tt,]) %*% yhat.Kalman[tt,]
      #logPL.DMA <-   sum(logPL*piMat)
    }
    
    #3f# DMS modelling
    thetahat.ms <- matrix(0, nrow = TT, ncol = ncol(X))
    logPL.DM.ms <- rep(0,TT)
    yhat.ms <- rep(0,TT)
    for (tt in 1:TT) {
      pimax <- which.max(piMat[tt,])
      thetahat.ms[tt,] <- thetahat[pimax,]
      logPL.DM.ms[tt] <- logPL[tt,pimax]
      #logPL.DMS <- sum(logPL.DM.ms)
    }  
    yhat.ms <- rowSums(X * thetahat.ms)
    
    #3g# MSFE
    Yt <- matrix(Y, ncol=1)
    MSE.DMA <- mean((Yt[(insample+1):TT,1] -  yhat.ma[(insample+1):TT])^2, na.rm = T)
    MSE.DMS <- mean((Yt[(insample+1):TT,1] -  yhat.ms[(insample+1):TT])^2, na.rm = T)
    MSE.OLS <- mean((Yt[(insample+1):TT,1] -  linear$fitted.values[(insample+1):TT])^2, na.rm = T)
    MSE.AR2 <- mean((AR$resid)^2, na.rm = T)
    
    #3h# MAFE 
    MAE.DMA <- mean(abs(Yt[(insample+1):TT,1] -  yhat.ma[(insample+1):TT]), na.rm = T)
    MAE.DMS <- mean(abs(Yt[(insample+1):TT,1] -  yhat.ms[(insample+1):TT]), na.rm = T)
    MAE.OLS <- mean(abs(Yt[(insample+1):TT,1] -  linear$fitted.values[(insample+1):TT]), na.rm = T)
    MAE.AR2 <- mean(abs(AR$resid), na.rm = T)
   
    #3i# logPL 
    logPL.DMA <- sum(logPL[(insample+1):TT,]*piMat[(insample+1):TT,])
    logPL.DMS <- sum(logPL.DM.ms[(insample+1):TT])
    
    
    list(MSE.DMA = MSE.DMA, MSE.DMS = MSE.DMS, MSE.OLS = MSE.OLS, MSE.AR2 = MSE.AR2, 
         MAE.DMA = MAE.DMA, MAE.DMS = MAE.DMS, MAE.OLS = MAE.OLS, MAE.AR2 = MAE.AR2, 
         logPL.DMA = logPL.DMA,logPL.DMS = logPL.DMS)
  }
  
  
  # DMA-DMS MODELLING FOR EACH CASE OF MACROECONOMIC SERIES
  testG = vector("list", ncol(macro.google))
  PerformanceG <- matrix(0, nrow = ncol(macro.google), ncol = 10)
  rownames(PerformanceG) <- c("y = Export.G","y = Import.G","y = CPI.G", "y = Spread10Y-2Y.G","y = Unemployment.G", "y = IPI.G", "y = VNSE.G")
  for (i in 1:ncol(macro.google)){
    testG[[i]] <- dmG(Y=dependent[[i]],x=explanatoryX[[i]],z=explanatoryZ[[i]], 
                      lambda = lambda, gamma = gamma, kappa=kappa, insample = insample)
    PerformanceG[i,] <- round(unlist(testG[[i]]), digit = 5)
  }
  colnames(PerformanceG) = c("MSE.DMA", "MSE.DMS","MSE.OLS","MSE.AR2", "MAE.DMA", "MAE.DMS", "MAE.OLS","MAE.AR2", "logPL.DMA", "logPL.DMS")
  return(PerformanceG)
  return(ProbG)
}

  ##3B## With Google as Model Probabilities
DMGProb <- function(lambda,gamma,kappa,insample){
  dmGProb<- function(Y,X,Z, insample, lambda , gamma, kappa, omega){
    #3a# Linear regression
    linear <- lm(Y ~ X-1)
    m = ncol(X)
    variance.y<- summary(linear)$sigma
    
    #3b# Auto regression order 2 
    AR <- ar(Y, order.max = 2, method = "ols", demean = T)
    
    #3c# Time-Varying Parameters model 
    
    #3c.1# Constructing the matrix to indicate which variables enter the model
    #For X
    models.which <- matrix(unlist(unique(combn(rep(c(0,1),m),m,simplify=FALSE))),(2^m),m,byrow=TRUE)
    models.which <- models.which[-(which(rowSums(models.which)==0)),]
    #For Z
    models.whichZ <- models.which[,-(1:2)]
    
    #3c.2# Initial matrices
    TT <- length(Y)
    K <- nrow (models.which)
    yhat <- matrix (rep(0,TT*K), ncol=K)
    theta<- vector ("list",K)
    sigma<- vector ("list",K)
    R <- vector ("list",K)
    theta0 <- rep(0,ncol(X))
    
    eAll <- matrix (rep(0,TT*K), ncol=K)
    VAll <- rep(0,K)
    thetahat <-matrix (rep(0,K*(ncol(X))),nrow=K)  
    Vtheta <- matrix (rep(0,K*(ncol(X))),nrow=K)  
    yhat.Kalman <- matrix (rep(0,TT*K), ncol=K)
    predvar.Kalman <- rep(0,K)
    
    model.dim <- rep(0,K)    #Number of predictors in model 
    model.dimZ <- rep(0,K)   #Number of predictors in model 
    for (k in 1:K){
      model.dim[k] <- sum(as.numeric(models.which[k,]))
      model.dimZ[k] <- sum(as.numeric(models.whichZ[k,]))
    } 
    
    s2y <- variance.y # variance of Y
    s2X <- rep(0,(ncol(X)))
    for (a in 1:ncol(X)){
      s2X[a] <- var(X[,a])
    }
    for (k in 1:K){
      theta[[k]] <-linear$coefficients[as.logical( models.which[k,])]
      if (model.dim[k] == 1) 
        sigma[[k]] <- s2y/s2X[as.logical(models.which[k,])]
      else
        sigma[[k]] <- diag(s2y/s2X[as.logical(models.which[k,])])
    }
    
    
    #3c.3# TVP model
    rm.Kalman <-function (theta, sigma, V,lambda, X, Y, TT) {
      R = list()
      theta = list(); theta[[1]] <- theta0
      sigma = list(); sigma[[1]] <- sigma0
      update = list()
      e <- rep(0,TT)
      V <-  (rep(0,TT));V[1]<- V0
      A <-  (rep(0,TT))
      yhat <- rep(0,TT)
      for (tt in 2:(TT)){
        R[[tt]] <- sigma[[tt-1]]/lambda
        if(ncol(as.matrix(sigma0))==1){
          xtmp <- X[tt]
        } else{
          xtmp <- X[tt,]
        }
        update[[tt]] <- R[[tt]]%*% xtmp %*%  solve( V[tt] + t(xtmp) %*% R[[tt]] %*% xtmp)
        yhat[tt] <- t(xtmp) %*% theta[[tt-1]]
        predvar <- t(xtmp) %*% R[[tt]] %*% xtmp + V[tt]
        e[tt] <- Y[tt]-yhat[tt]
        theta[[tt]] <- theta[[tt-1]] +update[[tt]] %*% e[tt]
        sigma[[tt]] <- R[[tt]] - update[[tt]] %*% t(xtmp) %*% R[[tt]]
        V[tt+1] <- kappa*V[tt] + (1-kappa)*e[tt]^2
        ifelse(V[tt+1] <= 0,V[tt],V[tt+1])
      }
      list (theta=theta[[TT]], sigma=sigma[[TT]],e=e, yhat=yhat, predvar=predvar, V = V[TT],R=R[[TT]])
    }
    
    #3c.3# K-cases of TVP model
    Kalmanall = vector("list",K)
    for (k in 1:K){
      if (model.dim[k] == 1){
        sigma0 <- s2y/s2X[as.logical( models.which[k,])]
      } else {
        sigma0 <- diag(s2y/s2X[as.logical( models.which[k,])]) 
      }
      theta0 <- linear$coefficients[as.logical( models.which[k,])]
      V0 <- variance.y<- summary(linear)$sigma
      Xtmp <-  X[,as.logical( models.which[k,])]
      Kalmanall[[k]] <-rm.Kalman (theta[[k]], sigma[[k]], V,lambda, Xtmp, Y, TT)
      ###
      thetahat[k,as.logical( models.which[k,])] <- Kalmanall[[k]]$theta
      Vtheta[k,as.logical( models.which[k,])] <- diag(Kalmanall[[k]]$sigma)
      yhat.Kalman[,k] <- Kalmanall[[k]]$yhat
      predvar.Kalman[k] <- Kalmanall[[k]]$predvar
      eAll[,k] <- Kalmanall[[k]]$e
      VAll[k] <- Kalmanall[[k]]$V
    } 
    
    
    #3d.A# Model Probabilities WITHOUT Google Probabilities
    eps=.001/nrow(models.which)
    model.update <- function (piold, gamma, eps,  y, yhat, predvar) {
      pipred <- (piold^gamma)+eps / sum((piold^gamma)+eps)
      logpyt <- -0.5*log(predvar) - 0.5*(y-yhat)^2/predvar # Update pi values
      pyt <- exp (logpyt)
      pinew <- pipred * pyt
      pinew <- pinew/sum(pinew)
      pinew <- pinew + eps
      pinew <- pinew/sum(pinew)
      
      list (pinew=as.vector(pinew), logpyt = logpyt)
    }
    
    pimat <- list()
    piOld <- rep (1/K, K)
    piMat <- matrix (rep(0,(TT)*K), ncol=K)
    logPL <- matrix (rep(0,(TT)*K), ncol=K)
    
    for (tt in 1:(TT)){
      # Update model probabilities
      pimat[[tt]] <- model.update(piOld, gamma, eps, Y[tt], yhat.Kalman[tt,],predvar.Kalman)
      piMat[tt,] <- pimat[[tt]]$pinew
      logPL[tt,] <- pimat[[tt]]$logpyt
    }
    
    
    #3d.II# Model Probabilities WITH Google Probabilities
    model.piGoogle <- function(i){
      if(model.dimZ[i] == 1)
        p1 = cbind(Z[,as.logical(models.whichZ[i,])],1)
      else
        p1 <-     Z[,as.logical(models.whichZ[i,])]
      if(model.dimZ[i] == 5)
        p2 = cbind(Z[,!as.logical(models.whichZ[i,])],1)
      else
        p2 <- 1 - Z[,!as.logical(models.whichZ[i,])]
      p1 <- apply(p1, 1, prod)
      p2 <- apply(p2, 1, prod)
      return(p1*p2)
    }
    
    piGoogle <- matrix(rep(0,(TT)*K), ncol=K)
    piGooglet <- matrix (rep(0,(TT)*K), ncol=K)
    for (k in (1:K)){
      # Update model probabilities with Google Prob.
      piGooglet[,k] <- unlist(lapply(k,model.piGoogle))
      piGoogle <- piGooglet/rowSums(piGooglet)
      piSum       <- omega*piMat + (1-omega)*piGoogle
    }
    
    #3e# DMA modelling
    thetahat.ma.piG <- matrix(0, nrow = TT, ncol = ncol(X))
    yhat.ma.piG <- rep(0,TT)
    for (tt in 1:(TT)){
      thetahat.ma.piG[tt,] <- t(piSum[tt,]) %*% thetahat
      yhat.ma.piG[tt] <- t(piSum[tt,]) %*% yhat.Kalman[tt,]
      #logPL.DMA.piG <-   sum(logPL*piSum)
    }
    
    #3f# DMS modelling
    thetahat.ms.piG <- matrix(0, nrow = TT, ncol = ncol(X))
    logPL.DM.ms.piG <- rep(0,TT)
    yhat.ms.piG <- rep(0,TT)
    for (tt in 1:TT) {
      pimax <- which.max(piSum[tt,])
      thetahat.ms.piG[tt,] <- thetahat[pimax,]
      logPL.DM.ms.piG[tt] <- logPL[tt,pimax]
      #logPL.DMS.piG <- sum(logPL.DM.ms.piG)
    } 
    yhat.ms.piG <- rowSums(X * thetahat.ms.piG)
    
    #3g# MSFE
    Yt <- matrix(Y, ncol=1)
    
    MSE.DMA.piG <- mean((Yt[(insample+1):TT,1] -  yhat.ma.piG[(insample+1):TT])^2, na.rm = T)
    MSE.DMS.piG <- mean((Yt[(insample+1):TT,1] -  yhat.ms.piG[(insample+1):TT])^2, na.rm = T)
    MSE.OLS <- mean((Yt[(insample+1):TT,1] -  linear$fitted.values[(insample+1):TT])^2, na.rm = T)
    MSE.AR2 <- mean((AR$resid)^2, na.rm = T)
    
    
    #3h# MAFE 
    MAE.DMA.piG <- mean(abs(Yt[(insample+1):TT,1] -  yhat.ma.piG[(insample+1):TT]), na.rm = T)
    MAE.DMS.piG<- mean(abs(Yt[(insample+1):TT,1] -  yhat.ms.piG[(insample+1):TT]), na.rm = T)
    MAE.OLS <- mean(abs(Yt[(insample+1):TT,1] -  linear$fitted.values[(insample+1):TT]), na.rm = T)
    MAE.AR2 <- mean(abs(AR$resid), na.rm = T)

    #3i# logPL 
    logPL.DMA.piG <- sum(logPL[(insample+1):TT,]*piSum[(insample+1):TT,])
    logPL.DMS.piG <- sum(logPL.DM.ms.piG[(insample+1):TT])
    
    list(MSE.DMA.piG=MSE.DMA.piG, MSE.DMS.piG=MSE.DMS.piG, MSE.OLS = MSE.OLS, MSE.AR2 = MSE.AR2,
         MAE.DMA.piG=MAE.DMA.piG, MAE.DMS.piG=MAE.DMS.piG, MAE.OLS = MAE.OLS, MAE.AR2 = MAE.AR2,  
         logPL.DMA.piG=logPL.DMA.piG, logPL.DMS.piG=logPL.DMS.piG )
  }
  
  # DMA-DMS MODELLING FOR EACH CASE OF MACROECONOMIC SERIES
  testpiG_1 = vector("list", ncol(macro.google))
  final.piG_1 <- matrix(0, nrow = ncol(macro.google), ncol = 10)
  
  testpiG_05 = vector("list", ncol(macro.google))
  final.piG_05 <- matrix(0, nrow = ncol(macro.google), ncol = 10)
  
  testpiG_00 = vector("list", ncol(macro.google))
  final.piG_00 <- matrix(0, nrow = ncol(macro.google), ncol = 10)
  
  rownames(final.piG_1)  = c("y = Export_1","y = Import_1","y = CPI_1", "y = Spread10Y-2Y_1","y = Unemployment_1", "y = IPI_1", "y = VNSE_1")
  rownames(final.piG_05) = c("y = Export_05","y = Import_05","y = CPI_05", "y = Spread10Y-2Y_05","y = Unemployment_05", "y = IPI_05", "y = VNSE_05")
  rownames(final.piG_00) = c("y = Export_00","y = Import_00","y = CPI_00", "y = Spread10Y-2Y_00","y = Unemployment_00", "y = IPI_00", "y = VNSE_00")
  
  for (i in 1:ncol(macro.google)){
    testpiG_1[[i]] <- dmGProb(Y=dependent[[i]], X=explanatoryX[[i]], Z=google.prob[[i]], 
                              lambda = lambda, gamma = gamma, kappa=kappa, insample = insample, omega = 1)
    final.piG_1[i,] <- round(unlist(testpiG_1[[i]]), digit = 5)
  }
  for (i in 1:ncol(macro.google)){
    testpiG_05[[i]] <- dmGProb(Y=dependent[[i]], X=explanatoryX[[i]], Z=google.prob[[i]], 
                               lambda = lambda, gamma = gamma, kappa=kappa, insample = insample, omega = 0.5)
    final.piG_05[i,] <- round(unlist(testpiG_05[[i]]), digit = 5)
  }
  for (i in 1:ncol(macro.google)){
    testpiG_00[[i]] <- dmGProb(Y=dependent[[i]], X=explanatoryX[[i]], Z=google.prob[[i]], 
                               lambda = lambda, gamma = gamma, kappa=kappa, insample = insample, omega = 0)
    final.piG_00[i,] <- round(unlist(testpiG_00[[i]]), digit = 5)
  } 
  all = interleave(final.piG_1,final.piG_05,final.piG_00)
  colnames(all) = c("MSE.DMA", "MSE.DMS","MSE.OLS","MSE.AR2", "MAE.DMA", "MAE.DMS", "MAE.OLS","MAE.AR2", "logPL.DMA", "logPL.DMS")
  return(all) 
}
    #Merging outcomes row-by-row
install.packages("gdata")
library(gdata)

##4## NOWCASTING PERFORMANCE
#Merging outcomes

ffactor95      <- DM(lambda=0.95,gamma=0.95,kappa=0.97,insample=48)
ffactor95GProb <- DMGProb(lambda=0.95,gamma=0.95,kappa=0.97,insample=48)

ffactor99      <- DM(lambda=0.99,gamma=0.99,kappa=0.97,insample=48)
ffactor99GProb <- DMGProb(lambda=0.99,gamma=0.99,kappa=0.97,insample=48)


##Export data##

#Gvar <- interleave(ffactor95,ffactor99)
#Gprob <- interleave(ffactor95GProb,ffactor99GProb)
#colnames(Gvar) = c("MSE.DMA", "MSE.DMS","MSE.OLS","MSE.AR2", "MAE.DMA", "MAE.DMS", "MAE.OLS","MAE.AR2", "logPL.DMA", "logPL.DMS")
#colnames(Gprob) = c("MSE.DMA", "MSE.DMS","MSE.OLS","MSE.AR2", "MAE.DMA", "MAE.DMS", "MAE.OLS","MAE.AR2", "logPL.DMA", "logPL.DMS")

#write.csv(Gvar,"~/Desktop/Gvar48.csv")
#write.csv(Gprob,"~/Desktop/Gprob48.csv")




####################
##### PLOTTING #####
####################
DMplot <- function(lambda,gamma,kappa,insample){
  dmG<- function(Y, x, z, insample, lambda , gamma, kappa){
    #3a# Linear regression
    linear <- lm(Y ~ x+z-1)
    m = ncol(x)
    variance.y<- summary(linear)$sigma
    
    #3b# Auto regression order 2 
    AR <- ar(Y, order.max = 2, method = "ols", demean = T)
    
    #3c# Time-Varying Parameters model 
    
    #3c.1Constructing the matrix to indicate which variables enter the model
    #For X
    models.whichX <- matrix(unlist(unique(combn(rep(c(0,1),m),m,simplify=FALSE))),(2^m),m,byrow=TRUE)
    models.whichX <- models.whichX[-(which(rowSums(models.whichX)==0)),]
    #For Z: Google variables conresponding to Macro. variables which enter the model 
    models.whichZ <- models.whichX[,-(1:2)]
    models.which <- cbind(models.whichX,models.whichZ)
    
    #3c.2# Initial matrices
    TT <- length(Y)
    K <- nrow (models.which)
    yhat <- matrix (rep(0,TT*K), ncol=K)
    theta<- vector ("list",K)
    sigma<- vector ("list",K)
    R <- vector ("list",K)
    theta0 <- rep(0,(ncol(x)+ncol(z)))
    
    eAll <- matrix (rep(0,TT*K), ncol=K)
    VAll <- rep(0,K)
    thetahat <-matrix (rep(0,K*(ncol(x)+ncol(z))),nrow=K)  
    Vtheta <- matrix (rep(0,K*(ncol(x)+ncol(z))),nrow=K)  
    yhat.Kalman <- matrix (rep(0,TT*K), ncol=K)
    predvar.Kalman <- rep(0,K)
    
    model.dim = rmAR <- rep(0,K)
    
    for (k in 1:K){
      model.dim[k] <- sum(as.numeric(models.which[k,]))
      rmAR[k] <- sum(as.numeric(models.which[k,-(1:2)]))
    }
    
    
    s2y <- variance.y # variance of Y
    s2x <- rep(0,(ncol(x)))
    s2z <- rep(0,(ncol(z)))
    for (a in 1:ncol(x)){
      for (b in 1:ncol(z)){
        s2x[a] <- var(x[,a])
        s2z[b] <- var(z[,b])
      }
    }
    s2X = c(s2x,s2z)
    for (k in 1:K){
      theta[[k]] <-linear$coefficients[as.logical( models.which[k,])]
      if (model.dim[k] == 1) 
        sigma[[k]] <- s2y/s2X[as.logical(models.which[k,])]
      else
        sigma[[k]] <- diag(s2y/s2X[as.logical(models.which[k,])])
    }
    #3c.3# TVP model
    rm.Kalman <-function (theta, sigma, V,lambda, X, Y, TT) {
      R = list()
      theta = list(); theta[[1]] <- theta0
      sigma = list(); sigma[[1]] <- sigma0
      update = list()
      e <- rep(0,TT)
      V <-  (rep(0,TT));V[1]<- V0
      A <-  (rep(0,TT))
      yhat <- rep(0,TT)
      for (tt in 2:(TT)){
        R[[tt]] <- sigma[[tt-1]]/lambda
        if(ncol(as.matrix(sigma0))==1){
          xtmp <- X[tt]
        } else{
          xtmp <- X[tt,]
        }
        update[[tt]] <- R[[tt]]%*% xtmp %*%  solve( V[tt] + t(xtmp) %*% R[[tt]] %*% xtmp)
        yhat[tt] <- t(xtmp) %*% theta[[tt-1]]
        predvar <- t(xtmp) %*% R[[tt]] %*% xtmp + V[tt]
        e[tt] <- Y[tt]-yhat[tt]
        theta[[tt]] <- theta[[tt-1]] +update[[tt]] %*% e[tt]
        sigma[[tt]] <- R[[tt]] - update[[tt]] %*% t(xtmp) %*% R[[tt]]
        V[tt+1] <- kappa*V[tt] + (1-kappa)*e[tt]^2
        ifelse(V[tt+1] <= 0,V[tt],V[tt+1])
      }
      list (theta=theta[[TT]], sigma=sigma[[TT]],e=e, yhat=yhat, predvar=predvar, V = V[TT],R=R[[TT]])
    }
    
    #3c.3# K-cases of TVP model
    Kalmanall = vector("list",K)
    X = cbind(x,z)
    for (k in 1:K){
      if (model.dim[k] == 1){
        sigma0 <- s2y/s2X[as.logical( models.which[k,])]
      } else {
        sigma0 <- diag(s2y/s2X[as.logical( models.which[k,])]) 
      }
      theta0 <- linear$coefficients[as.logical( models.which[k,])]
      V0 <- variance.y<- summary(linear)$sigma
      Xtmp <-  X[,as.logical( models.which[k,])]
      Kalmanall[[k]] <-rm.Kalman (theta[[k]], sigma[[k]], V,lambda, Xtmp, Y, TT)
      ###
      thetahat[k,as.logical( models.which[k,])] <- Kalmanall[[k]]$theta
      Vtheta[k,as.logical( models.which[k,])] <- diag(Kalmanall[[k]]$sigma)
      yhat.Kalman[,k] <- Kalmanall[[k]]$yhat
      predvar.Kalman[k] <- Kalmanall[[k]]$predvar
      eAll[,k] <- Kalmanall[[k]]$e
      VAll[k] <- Kalmanall[[k]]$V
    } 
    
    #3d# Model Probabilities
    eps=.001/nrow(models.which)
    model.update <- function (piold, gamma, eps,  y, yhat, predvar) {
      pipred <- (piold^gamma)+eps / sum((piold^gamma)+eps)
      logpyt <- -0.5*log(predvar) - 0.5*(y-yhat)^2/predvar # Update pi values
      pyt <- exp (logpyt)
      pinew <- pipred * pyt
      pinew <- pinew/sum(pinew)
      pinew <- pinew + eps
      pinew <- pinew/sum(pinew)
      
      list (pinew=as.vector(pinew), logpyt = logpyt)
    }
    
    pimat <- list()
    piOld <- rep (1/K, K)
    piMat <- matrix (rep(0,(TT)*K), ncol=K)
    logPL <- matrix (rep(0,(TT)*K), ncol=K)
    #3d.1# Update model probabilities
    for (tt in 1:(TT)){
      pimat[[tt]] <- model.update(piOld, gamma, eps, Y[tt], yhat.Kalman[tt,],predvar.Kalman)
      piMat[tt,] <- pimat[[tt]]$pinew
      logPL[tt,] <- pimat[[tt]]$logpyt
    }
    Ssize = rep(0,TT)
    Prob = matrix(0,ncol=12, nrow=TT)
    for (tt in 1:TT) {
      pimax <- which.max(piMat[tt,])
      Ssize[tt] <- rmAR[pimax]
      #logPL.DMS <- sum(logPL.DM.ms)
      Prob[tt,] <- piMat[tt,pimax]*models.which[pimax,-(1:2)]
      
    }  
    #Prob <- piMat%*%models.which[,-(1:2)]
    Esize <- piMat%*%rmAR
    
    list(Prob = Prob, Esize = Esize, Ssize=Ssize)
  }
  
  
  # DMA-DMS MODELLING FOR EACH CASE OF MACROECONOMIC SERIES
  testG = vector("list", ncol(macro.google))
  
  #rownames(PerformanceG) <- c("y = Export.G","y = Import.G","y = CPI.G", "y = Spread10Y-2Y.G","y = Unemployment.G", "y = IPI.G", "y = VNSE.G")
  for (i in 1:ncol(macro.google)){
    testG[[i]] <- dmG(Y=dependent[[i]],x=explanatoryX[[i]],z=explanatoryZ[[i]], 
                      lambda = lambda, gamma = gamma, kappa=kappa, insample = insample)
  }
  #colnames(PerformanceG) = c("MSE.DMA", "MSE.DMS","MSE.OLS","MSE.AR2", "MAE.DMA", "MAE.DMS", "MAE.OLS","MAE.AR2", "logPL.DMA", "logPL.DMS")
  return(testG)
}

ProbG      <- DMplot(lambda=0.99,gamma=0.99,kappa=0.97,insample=48)

##2## Number of Predictors 

graphics.off()
par(mar=c(2,2,2,2))
par(mfrow=c(7,1))

Esize.all = Ssize.all <- matrix(0, ncol = ncol(macro.google), nrow = nrow(macro.google)-2 )

for (i in 1:ncol(macro.google)){
  Esize.all[,i] <-round(ProbG[[i]]$Esize,digit=3)
  Esize.all <- ts(Esize.all,end=end(macro.google),frequency = 12)
  Ssize.all[,i] <-round(ProbG[[i]]$Ssize,digit=3)
  Ssize.all <- ts(Ssize.all,end=end(macro.google),frequency = 12)
}
legend <- c("DMA","DMS")
main <- c("Export","Import","CPI", "Spread10Y-2Y","Unemployment", "IPI", "VNSE")

for (i in 1:ncol(macro.google)){
  #jpeg(file="~/Desktop/PLOT/7.jpeg", width=800, height=400)
  plot(Esize.all[,i], col= 1,ylim = c(1,12),type="l",lwd=3, ylab = "Number of Predictors", xlab = "Date")
  lines(Ssize.all[,i],col= 1)
  axis(2, at = 1:12, tck = 1, lty = 1,lwd=.1, col = "grey", labels = NA)
  legend("topright",legend, col=1, lwd = c(3,1), cex=.9)
  #dev.off()
}

##2## Posterior probabilities

graphics.off()
par(mar=c(2,2,2,2))
par(mfrow=c(4,2))

library(ggplot2)
Prob  <- matrix(0, nrow = nrow(macro.google)-2, ncol = 12)
cname <- c("Export", "Import","CPI", "Spread10Y-2Y","Unemployment", "IPI", "VNSE",
           "ExportG","Import.G","CPI.G", "Spread10Y-2Y.G","Unemployment.G", "IPI.G", "VNSE.G")
Date <- seq(as.POSIXct("2004-03-01"),as.POSIXct("2019-12-01"),by="1 month")

for (i in 1:ncol(macro.google)){
  Prob <-round(ProbG[[i]]$Prob,digit=3)
  colnames(Prob) <- cname[-c(i,i+7)]
  df <- data.frame(Prob)
  df[df == 0] <- NA
  df_with_time <- data.frame(Date,df)
  data <- reshape2::melt(df_with_time,id.vars="Date",variable.name = 'Variables',value.name = 'Prob.')
  #pdf(file="~/Desktop/PLOT/P_7.pdf")
  ggplot(data,aes(x = Date,y = Variables, size = Prob.)) +geom_line() +
    theme_bw()+
    theme(
      legend.title = element_text(color = "black", size = 10),
      legend.text = element_text(color = "black", size = 8),
      panel.background = element_rect(colour = "white"),
      panel.grid.major = element_line(colour = "black",size = 0.1)
    )
  #dev.off()
}  





