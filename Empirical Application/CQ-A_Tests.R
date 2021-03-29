## First Test


#' compute (gloabl) quantile test statistic (cf. Section 3)
#'
#' @param Y n-dim. vector containing the observations on the outcome
#' @param n
#' @param IndX 
#' @param IndZ n-dim. vector containing the observations on the estimated propensity score
#' @param qhat d-dim. vector containing the evluation quantiles
#' @param taueval 
VBDNStat <- function(Y,n, IndX, IndZ, qhat, taueval) {

  
  # nonparametric quantile regression of Y on X
  Tn.vals <- vector(mode="list",dim(as.matrix(taueval))[1])
  epsilonhat <- vector(mode="list",dim(as.matrix(taueval))[1])
  
  for (tau in 1:dim(as.matrix(taueval))[1]) {
    
    epsilonhat[[tau]] <- (Y - qhat[[tau]])
    
    # test statistic
      for (i in 1:dim(IndX)[2]) {
        Tn.vals[[tau]][i] <- mean(((epsilonhat[[tau]] <=0)-taueval[[tau]])*IndX[,i]*IndZ[,i]) 
      }
  
  }
  
  teststat <- max(abs(sqrt(n)*unlist(Tn.vals)))
  
  return(list(teststat=teststat, epsilonhat=epsilonhat))
}


#' compute boostrap statistic
#'
#' @param Bernl n-dim. vector containing simulated Bernoulli random variables
#' @param n 
#' @param IndX 
#' @param IndZ 
#' @param taueval 
#' @param stat
#' WBootStat(Bernl,n,IndX,IndZ,taueval,stat)
WBootStat <- function(Bernl,n, IndX, IndZ, taueval,stat) {
  
  
  Tn.vals  <- vector(mode="list",dim(as.matrix(taueval))[1])
  
  for (tau in 1:dim(as.matrix(taueval))[1]) {
    
      if (stat==1) {
        for (i in 1:dim(IndX)[2]) {
          Tn.vals[[tau]][i] <- mean((Bernl[[tau]] - taueval[tau])*IndX[,i]*IndZ[,i])
        }    
      } else {
        for (i in 1:dim(IndX)[2]) {
          Tn.vals[[tau]][i] <- mean((Bernl[[tau]] - taueval[tau])*IndX[,i]*IndZ[[tau]][,i])
        }
      }
  }
  
  
  teststatb <- max(abs(sqrt(n)*unlist(Tn.vals)))
  
  return(list(teststat=teststatb))
}



#' perform (local) quantile test statistic (cf. Section 4)
#'	
#' compute test statistic for testing the hypothesis H0: Pr(Pr(Y<= Q_tau[Y|X] |P]=tau) = 1
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param IndX 
#' @param IndZ 
#' @param T1IndZ
#' @param qhat 
#' @param taueval 
#' @param size 
#' @param B 
#' @param stat 
#' VBDNTest(Y, X, Z, delta=.75,hp=.25,size=0.05, taueval=c(.25,.5,.75), B=100, a=NA, ckertype="epanechnikov")
VBDNTest <- function(Y, IndX, IndZ, T1IndZ, qhat, taueval, size, B, stat) {
  
  stopifnot(size<1 & size>0)
  n <- length(Y)
  # compute test statistic
  Tstat <- VBDNStat(Y,n, IndX, IndZ, qhat, taueval)
  teststat <- Tstat$teststat
  epsilonhat <- Tstat$epsilonhat
  
  
  teststatb <- rep(0,B)
  Bernl <- vector(mode="list",dim(as.matrix(taueval))[1])
  for (b in 1:B) {
    
  Bernoullirv = runif(length(Y), min = 0, max = 1)
  
  for (tau in 1:dim(as.matrix(taueval))[1]) {  
  Bernl[[tau]] = (Bernoullirv <= taueval[[tau]])*1
  }
  # bootstrap critical value
  teststatb[b] <- WBootStat(Bernl,n , IndX, T1IndZ, taueval, stat)$teststat

  }
  cv <- quantile(teststatb, 1-size, na.rm=TRUE)
  
  Fn <- ecdf(teststatb)
  pval <- 1-Fn(teststat)
  
  return(list(teststat=teststat, cv=cv, rej=teststat>cv, pval=pval))
}

## Second Test


#' compute localized Volgushev et al. (2013) test statistic
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param qhat
#' @param kweights
#' @param taueval)
#' LocVBDNStat(qhat, kweights,taueval)
LocVBDNStat <- function(Y, qhat, kweights ,taueval) {


  # nonparametric quantile regression of Y on X
  Tn.vals <- vector(mode="list",dim(as.matrix(taueval))[1])
  Var.vals <- vector(mode="list",dim(as.matrix(taueval))[1])
  epsilonhat <- vector(mode="list",dim(as.matrix(taueval))[1])
  
  
  for (tau in 1:dim(as.matrix(taueval))[1]) {
    
    Tn.vals[[tau]] <- sum((((Y - qhat[[tau]]) <=0) - taueval[tau])*kweights)
    Var.vals[[tau]] <- (3/5)*sum(kweights*((((Y - qhat[[tau]]) <=0) - taueval[tau])^2))
    epsilonhat[[tau]] <- (Y - qhat[[tau]])
    
  }
  
  teststat <- max(abs((unlist(Tn.vals))/sqrt(unlist(Var.vals))),na.rm=TRUE)
  
  return(list(teststat=teststat, epsilonhat=epsilonhat))
}


#' compute boostrap statistic
#'
#' @param Bernl n-dim. vector containing simulated Bernoulli random variables
#' @param kweights
#' @param taueval
#' 
#' WBootStat(Bernl,kweights , taueval)
LocWBootStat <- function(Bernl,kweights , taueval) {
  
  Tn.vals <- vector(mode="list",dim(as.matrix(taueval))[1])
  Var.vals <- vector(mode="list",dim(as.matrix(taueval))[1])
  
  for (tau in 1:dim(as.matrix(taueval))[1]) {
    
    # bootstrap statistic
        Tn.vals[[tau]] <- sum((Bernl[[tau]] - taueval[tau])*kweights)
        Var.vals[[tau]] <- (mean((Bernl[[tau]] - taueval[tau])^2))*(sum(kweights))*(3/5)          
        
   } 
  # bootstrap statistic
  
  teststatb <- max(abs(unlist(Tn.vals)/sqrt(unlist(Var.vals))),na.rm=TRUE)
  
  return(list(teststat=teststatb))
}



#' perform localized VBDN test to discriminate between selection and misspecification
#'	
#' compute test statistic for testing the hypothesis H0: Pr(Pr(Y<= Q_tau[Y|X] |X,P=1]=tau) = 1
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param qhat
#' @param kweights
#' @param size
#' @param taueval
#' @param B
#' LocVBDNTest(Y, qhat, kweights ,size, taueval, B)
LocVBDNTest <- function(Y, qhat, kweights ,size, taueval, B) {

 
  
  stopifnot(size<1 & size>0)


  # compute test statistic
  Tstat <- LocVBDNStat(Y,qhat, kweights, taueval)
  teststat <- Tstat$teststat
  epsilonhat <- Tstat$epsilonhat

  
  teststatb <- rep(0,B)
  Bernl <- vector(mode="list",dim(as.matrix(taueval))[1])
  for (b in 1:B) {
    for (tau in 1:dim(as.matrix(taueval))[1]) {  
      Bernl[[tau]] = (runif(length(Y), min = 0, max = 1) <= taueval[[tau]])*1
    }
    # bootstrap critical value
    teststatb[b] <- LocWBootStat(Bernl, kweights, taueval)$teststat
    
  }
  cv <- quantile(teststatb, 1-size, na.rm=TRUE)
  
  Fn <- ecdf(teststatb)
  pval <- 1-Fn(teststat)
  
  return(list(teststat=teststat, cv=cv, rej=teststat>cv, pval=pval))
}


