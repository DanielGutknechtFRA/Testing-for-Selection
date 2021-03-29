## First Test

#' draw a random sample from Mammen's two-point distribution
#'	
#' @param n sample size
#' rMammen(n)
rMammen <- function(n) {
  V <- rbinom(n, 1, prob=(sqrt(5)-1)/(2*sqrt(5)))
  V[V==0] <- (1-sqrt(5))/2
  V[V==1] <- (1+sqrt(5))/2
  return(V)
}


#' compute Delgado and Manteiga (2001) test statistic
#'
#' compute Delgado and Manteiga (2001) test statistic for testing the hypothesis H0: E[Y|X,Z] = E[Y|X]
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param X matrix with n rows containing the observations on the scalar or vector X
#' @param Z matrix with n rows containing the observations on the scalar or vector Z
#' @param a vector of bandwidths, of the same dimension as there are columns in X, if unspecified, then the bandwidths are determined by cross-validation from nonparametric regression of Y on X
#' @param ckertype character string denoting the kernel function to be used, as in np package (default: "gaussian")
#' @param stat character string denoting the type of test statistic to be computed: Cramer-von-Mises ("CvM", default) or Kolmogorov-Smirnov ("KS")
#' computeDMStat(Y, X, Z, a=NA, ckertype="gaussian", stat="CvM")
DMStat <- function(Y, X, Z, a, ckertype, stat) {
  
  n <- length(Y)
  if (is.matrix(X) | is.data.frame(X)) dX <- ncol(X) else { stopifnot(is.numeric(X)); dX <- 1; }
  if (is.matrix(Z) | is.data.frame(Z)) dZ <- ncol(Z) else { stopifnot(is.numeric(Z)); dZ <- 1; }
  
  # nonparametric regression of Y on X
  fhat <- c(np::npksum(txdat=X, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / n
  Yhat <- c(np::npksum(txdat=X, tydat=Y, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / fhat / n
  epsilonhat <- Y-Yhat
  
  # test statistic
  if (dX==1 & dZ==1) {
    Tn.vals <- rep(NA,n)
    for (i in 1:n) {
      ind <- (X<=X[i])*(Z<=Z[i])
      Tn.vals[i] <- mean(epsilonhat*fhat*ind)
    }
  } else {
    Tn <- function(x,z) {
      X <- as.matrix(X); Z <- as.matrix(Z)
      ind <- apply((X<=x), 1, prod) * apply((Z<=z), 1, prod)
      return(mean(epsilonhat*fhat*ind))
    }
    Tn.vals <- c(apply(cbind(X,Z), 1, function(x) Tn(x[1:dX],x[(dX+1):(dX+dZ)])))
  }
  
  teststat <- switch(stat,
                     "CvM" = sum(Tn.vals^2),
                     "KS"  = max(abs(sqrt(n)*Tn.vals)))
  
  return(list(teststat=teststat, epsilonhat=epsilonhat, Yhat=Yhat))
}


#' perform the Delgado and Manteiga (2001) test
#'
#' perform the Delgado and Manteiga (2001) test of the hypothesis H0: E[Y|X,Z] = E[Y|X]
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param X matrix with n rows containing the observations on the scalar or vector X
#' @param Z matrix with n rows containing the observations on the scalar or vector Z
#' @param size scalar between 0 and 1, denoting the nominal size of the test (default: 0.05)
#' @param B integer denoting the number of bootstrap samples to be used (default: 100)
#' @param a vector of bandwidths, of the same dimension as there are columns in X, if unspecified, then the bandwidths are determined by cross-validation from nonparametric regression of Y on X
#' @param ckertype character string denoting the kernel function to be used, as in np package (default: "gaussian")
#' @param stat character string denoting the type of test statistic to be computed: Cramer-von-Mises ("CvM", default) or Kolmogorov-Smirnov ("KS")
#' DMTest(Y, X, Z, size=0.05, B=100, a=NA, ckertype="gaussian", stat="KS")
DMTest <- function(Y, X, Z, size, B, a, ckertype, stat) {
  
  stopifnot(size<1 & size>0)
  
  # compute test statistic
  res <- DMStat(Y, X, Z, a, ckertype, stat)
  teststat <- res$teststat
  epsilonhat <- res$epsilonhat
  Yhat <- res$Yhat

  # bootstrap critical value
  teststatb <- rep(0,B)
  for (b in 1:B) {
    V <- rMammen(length(Y))
    Yhatstar <- Yhat+epsilonhat*V
    teststatb[b] <- DMStat(Yhatstar, X, Z, a, ckertype, stat)$teststat
  }
  cv <- quantile(teststatb, 1-size, na.rm=TRUE)
  
  Fn <- ecdf(teststatb)
  pval <- 1-Fn(teststat)
  
  return(list(teststat=teststat, cv=cv, rej=teststat>cv, pval=pval, a=ahat))
}



## SECOND TEST

#' compute localized Delgado and Manteiga (2001) test statistic for testing the hypothesis H0: E[Y|X,Z] = E[Y|X]
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param X
#' @param Z
#' @param kweights
#' @param a
#' @param ckertype
#' @param stat

LocDMStat <- function(Y,X,Z,kweights,a,ckertype, stat) {
  
  n <- length(Y)
  if (is.matrix(X) | is.data.frame(X)) dX <- ncol(X) else { stopifnot(is.numeric(X)); dX <- 1; }
  if (is.matrix(Z) | is.data.frame(Z)) dZ <- ncol(Z) else { stopifnot(is.numeric(Z)); dZ <- 1; }
  
  # nonparametric regression of Y on X
  if (any(is.na(a)) | length(a)!=dX) a <- (np::npregbw(xdat=X, ydat=Y))$bw
  fhat <- c(np::npksum(txdat=X, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / n
  Yhat <- c(np::npksum(txdat=X, tydat=Y, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / fhat / n
  epsilonhat <- Y-Yhat

  # test statistic
  if (dX==1 & dZ==1) {
    Tn.vals <- rep(NA,n)
    for (i in 1:n) {
      ind <- (X<=X[i])*kweights
      Tn.vals[i] <- mean(epsilonhat*fhat*ind)
    }
  } else {
    Tn <- function(x,z) {
      X <- as.matrix(X); Z <- as.matrix(Z)
      ind <- apply((X<=x), 1, prod)*kweights
      return(mean(epsilonhat*fhat*ind))
    }
    Tn.vals <- c(apply(cbind(X,Z), 1, function(x) Tn(x[1:dX],x[(dX+1):(dX+dZ)])))
  }
  
  teststat <- switch(stat,
                     "CvM" = sum(Tn.vals^2),
                     "KS"  = max(abs(sqrt(n)*Tn.vals)))
  
  return(list(teststat=teststat, epsilonhat=epsilonhat, Yhat=Yhat))
}


#' perform local version of Delgado and Manteiga (2001) test
#'
#' perform local version of Delgado and Manteiga (2001) test of the hypothesis H0: E[Y|X,Z=1] = E[Y|X]
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param X
#' @param Z
#' @param delta
#' @param hp
#' @param size
#####' @param B####
#' @param a
#' @param ckertype
#' @param stat
#LocDMTest <- function(Y, X, Z, delta, hp, size, B, a, ckertype, stat) {
LocDMTest <- function(Y, X, Z, delta, hp, size, a, ckertype, stat) {
    
    
  stopifnot(size<1 & size>0)
  
  # compute kernel weights
  kweights <-  ((np::npksum(txdat=Z,exdat=delta, leave.one.out=FALSE, bws=hp, ckertype=ckertype, return.kernel.weights=TRUE)$kw)/hp)
  
  
  # compute test statistic
  res <- LocDMStat(Y, X, Z, kweights,a,ckertype, stat)
  teststat <- res$teststat
  epsilonhat <- res$epsilonhat
  Yhat <- res$Yhat

  # bootstrap critical value
  #teststatb <- rep(0,B)
  #for (b in 1:B) {
  #  V <- rMammen(length(Y))
  #  Yhatstar <- Yhat+epsilonhat*V
  #  teststatb[b] <- LocDMStat(Yhatstar, X, Z, kweights, a, ckertype, stat)$teststat
  #}
  #cv <- quantile(teststatb, 1-size, na.rm=TRUE)
  #Fn <- ecdf(teststatb)
  #pval <- 1-Fn(teststat)
  
  
  cv <- qnorm(1-(size/2))
  pval <- 1-pnorm(teststat,mean=0,sd=1)
  
  return(list(teststat=teststat, cv=cv, rej=abs(teststat)>cv, pval=pval, a))
}

