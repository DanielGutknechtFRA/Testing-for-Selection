## First Test

#' draw a random sample from Mammen's two-point distribution
#' rMammen(n)
rMammen <- function(n) {
  V <- rbinom(n, 1, prob=(sqrt(5)-1)/(2*sqrt(5)))
  V[V==0] <- (1-sqrt(5))/2
  V[V==1] <- (1+sqrt(5))/2
  return(V)
}


#' compute Delgado and Manteiga (2001) test statistic
DMStat <- function(Y,X, IndX, IndZ, a, ckertype, stat) {
  
  n <- length(Y)
  if (is.matrix(X) | is.data.frame(X)) dX <- ncol(X) else { stopifnot(is.numeric(X)); dX <- 1; }

  # nonparametric regression of Y on X
  if (any(is.na(a)) | length(a)!=dX) a <- (np::npregbw(xdat=X, ydat=Y))$bw
  fhat <- c(np::npksum(txdat=X, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / n
  Yhat <- c(np::npksum(txdat=X, tydat=Y, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / fhat / n
  epsilonhat <- Y-Yhat
  
  # test statistic
    Tn.vals <- rep(NA,n)
    for (i in 1:n) {
      Tn.vals[i] <- mean(epsilonhat*fhat*IndX[,i]*IndZ[,i])
    }

  teststat <- switch(stat,
                     "CvM" = sum(Tn.vals^2),
                     "KS"  = max(abs(sqrt(n)*Tn.vals)))
  
  return(list(teststat=teststat, epsilonhat=epsilonhat, Yhat=Yhat, a=a))
}


#' perform the Delgado and Manteiga (2001) test
DMTest <- function(Y,X, IndX, IndZ, size, B, a, ckertype, stat) {
  
  stopifnot(size<1 & size>0)
  
  # compute test statistic
  res <- DMStat(Y,X, IndX, IndZ, a, ckertype, stat)
  teststat <- res$teststat
  epsilonhat <- res$epsilonhat
  Yhat <- res$Yhat
  ahat <- res$a
  
  # bootstrap critical value
  teststatb <- rep(0,B)
  for (b in 1:B) {
    V <- rMammen(length(Y))
    Yhatstar <- Yhat+epsilonhat*V
    teststatb[b] <- DMStat(Yhatstar,X, IndX, IndZ, ahat, ckertype, stat)$teststat

  }

return(list(teststat=teststat, teststatb=teststatb))
}


