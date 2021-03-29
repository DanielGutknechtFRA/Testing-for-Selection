
#' compute localized Delgado and Manteiga (2001) test statistic for testing the hypothesis H0: E[Y|X,Z] = E[Y|X]

LocDMStat <- function(Testvar,kweights) {
  
  
  
  # test statistic
      Tn.vals <- sum(Testvar*kweights)
      Var.vals <- (3/5)*(sum((Testvar^2)*kweights))

  teststat <-  abs((unlist(Tn.vals))/sqrt(unlist(Var.vals)))
  
  return(list(teststat=teststat))
}


#' perform local version of Delgado and Manteiga (2001) test
LocDMTest <- function(Testvar, kweights, size) {
  
  n <- length(Testvar)

  
  # compute test statistic
  teststat <- LocDMStat(Testvar, kweights)

  # bootstrap critical value
  # bootstrap critical value
  #teststatb <- rep(0,B)
  #for (b in 1:B) {
  #  V <- rMammen(length(Testvar))
  #  
  #      Tn.vals <- sum(V*Testvar*kweights)
  #      Var.vals <- mean(V^2)*(3/5)*sum((Testvar^2)*kweights)

  #  teststatb[b] <- abs((unlist(Tn.vals))/sqrt(unlist(Var.vals)))
  #}
  #cv <- quantile(teststatb, 1-size, na.rm=TRUE)
  #Fn <- ecdf(teststatb)
  #pval <- 1-Fn(teststat)
  cv <- qnorm(1-(size/2))
  
  
  #return(list(teststat=teststat, cv=cv, rej=teststat>cv, pval=pval, a))
  return(list(teststat=teststat, cv=cv, rej=teststat>cv))
}

