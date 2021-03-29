########################################################
## Simulation Routine - Conditional Quantile - Test 2 ##
##         -   No Data-Driven Bandwidth  -            ##
########################################################

rm(list = ls())


library(np)
library(tidyr)
library(quantreg)
library(MASS)
source("Multiplot.R")


## Second Test
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
#' WBootStat(Bernl,kweights, taueval)
#' 
LocWBootStat <- function(Bernl,kweights , taueval) {
  
  Tn.vals <- vector(mode="list",dim(as.matrix(taueval))[1])
  Var.vals <- vector(mode="list",dim(as.matrix(taueval))[1])
  
  for (tau in 1:dim(as.matrix(taueval))[1]) {
    
    # bootstrap statistic
    Tn.vals[[tau]] <- sum((Bernl[[tau]] - taueval[tau])*kweights)
    Var.vals[[tau]] <- mean((Bernl[[tau]] - taueval[tau])^2)*sum(kweights)*(3/5)          
  } 
  # bootstrap statistic
  
  teststatb <- max(abs(unlist(Tn.vals)/sqrt(unlist(Var.vals))),na.rm=TRUE)
  
  return(list(teststat=teststatb))
}


#' perform Volgushev et al. (2013) test 
#'	
#' compute test statistic for testing the hypothesis H0: Pr(Pr(Y<= Q_tau[Y|X] |X,P]=tau) = 1
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param X matrix with n rows containing the observations on the scalar or vector X
#' @param Z n-dim. vector containing the observations on the estimated propensity score
#' @param size scalar between 0 and 1, denoting the nominal size of the test
#' @param taueval d-dim. vector of evaluation quantiles
#' @param B integer denoting the number of bootstrap samples to be used (default: 100)
#' @param a vector of bandwidths, of the same dimension as there are columns in X, if unspecified, then the bandwidths are determined by cross-validation from nonparametric regression of Y on X
#' @param ckertype character string denoting the kernel function to be used, as in np package (default: "gaussian")
#' @param stat use marginal or conditional cdf to construct bootstrap statistics
#' @return a list containing the following elements: 'teststat' value of the test statistic, 'cv' bootstrap critical value, 'rej' a 1-0 indicator for whether the test rejects or not, 'pval' p-value, 'a' the bandwidth(s)
#' @keywords discrimination selection vs misspecification, local test
#' @export
#' @examples
#' VBDNTest(Y, X, Z, delta=.75,hp=.25,size=0.05, taueval=c(.25,.5,.75), B=100, a=NA, ckertype="epanechnikov")

LocVBDNTest <- function(Y, qhat, kweights, taueval, B) {
    
  
  n <- length(Y)
  # compute test statistic
  
  
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
  
  return(list(teststat=teststat, teststatb=teststatb))
}


time=proc.time()

# Set working directory

setwd("C:/Users/gutknech/Downloads/Simulations/Results")

# Set number of Monte Carlo replications
MC <- 999

# Set the nominal size of both tests
alpha1 <- 0.05
alpha2 <- 0.1


CQ <- 1

# Set the number of bootstrap replications
Bsize <- 1


  deltavec <- c(0.95,0.95,0.975,.98)
  hpvec    <- c(0.075,0.05,0.03,.025)
  # Rejection Rate list
  RejRate <- lapply(1:length(hpvec), function(x) matrix(0,length(BW),4))
  


# Define the grid of quantiles
tau_grid = c(.3,.5,.7) 

# Use the marginal CDF of Phat (MD==1: `yes'; MD==0: `no')
MD <- 1

# Set the degree of selection 
rho <- 0.25

# Set the misspecification parameter:
gamma1 <- 0.25
gamma2 <- 1

T1 <- 1

# Set trimming constant
trim <- .025


# Selectoin Equation:
# 1 - continuous N(0,1) 
# 2 - bivariate Binomial {-0.5,0.5}
# 3 - discrete Poisson(1.5)
# 4 - discrete Uniform {-3,-2,-1,0,1,2,3}
#SE <- 2

# Set the sample sizes

SS <- matrix(c(1000, 2000),1,2)

# Results matrix
  BW<-c(3.5,4.0,4.5)
#  BW <- c(1)
  TestRes <- lapply(1:length(hpvec), function(x) matrix(0,MC,dim(SS)[2]*3))
  BootRes <- lapply(1:length(hpvec), function(x) matrix(0,MC,dim(SS)[2]*3))
  

## Evaluation matrix functions

IndMat <- function(Z) {
  Z <- as.matrix(Z)
  indz <- matrix(NA,dim(Z)[1],dim(Z)[1])
  for (i in 1:dim(Z)[1]) {
    indz[,i] <- as.numeric(Z<=Z[i])
  }
  return(indz)
}

T1IndMat <- function(Z,Fz) {
  Z <- as.matrix(Z)
  indz <- matrix(NA,dim(Z)[1],dim(Z)[1])
  
  for (i in 1:dim(Z)[1]) {
    indz[,i] <- ((Z<=Z[i])-Fz[i])
  }
  return(indz)
}



Kweights <- function(arg) {
  
  arg <- as.matrix(arg)
  return((3/4)*(1-(arg^2))*(as.numeric(abs(arg)<=1)))
  
}

mnpCV <- function(BWgrid,Y,X) {
  n <- dim(as.matrix(X))[1]
  CV_err_h =rep(NA,length(BWgrid))
  for(j in 1:length(BWgrid)){
  h_using = (BWgrid[j]*apply(X,2,sd)*(n)^(-1/6))
  fZhat <- np::npksum(txdat=X,leave.one.out=TRUE,  bws=h_using, ckertype="epanechnikov")$ksum
  Fhat <- np::npksum(txdat=X, tydat=as.numeric(Y), leave.one.out=TRUE, bws=h_using, ckertype="epanechnikov")$ksum / (fZhat+0.00001) 
  CV_err_h[j] = mean((Y - as.matrix(Fhat)^2))
}

return(BWgrid[which(CV_err_h ==min(CV_err_h))])

}
##########################
##########################
### SIMULATION ROUTINE ###
##########################
##########################

for (SE in 1:1) {
  
  # Set random number seed and get start time
  
  set.seed(12345)
  
  
  for (ss in 1:dim(SS)[2]) {
    
    
    Ssize <- SS[ss]
    
    for (mc in 1:MC) {
      ########################
      # Data Generation      #
      ########################
      
      if ((mc %in% seq(100,100*floor(MC/100),100))==TRUE){
        cat(mc)
      } else {
        cat('.')
      }
      
      
        
        if (SE==1) {
          Z <- rnorm(Ssize,mean=0, sd=1)
        } else if (SE==2) {
          Z <- 0.5-rbinom(Ssize, 1, 0.5)
        } 
        
        X <- runif(Ssize, min = 0, max = 1)
        E <- mvrnorm(Ssize,matrix(c(0,0),2,1),matrix(c(1,rho,rho,1),2,2,byrow=TRUE))
        S <- ((0.75*(X-0.5)+(0.75*Z))>=(gamma2*E[,2]))
        
        
        Y <- (X-0.5)^3 + (X-0.5)^2 + (X-0.5)+ (gamma1*Z)  + 0.5*E[,1]
          
        if (SE==1) {
          Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, (gamma2*1))
        } else if (SE==2) {
          Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, (gamma2*1))
        }    
        
      Y <- Y[S]
      X <- X[S]
      Pz <- Pz[S]
#      E1 <- E[S,1]
      
      filter <- (X<quantile(X,probs=trim))|(X>quantile(X,probs=1-trim))
      Xtrim <- X[!filter]
      Ptrim <- Pz[!filter]
      Ytrim <- Y[!filter]
#      Etrim <- E1[!filter]
      rm(filter)
      
      
      ########################
      # Conditional Quantile #
      ########################
        for (h in 1:length(BW)) {
          
          qhat <- vector(mode="list",dim(as.matrix(tau_grid))[1])
          
          for (tau in 1:dim(as.matrix(tau_grid))[1]) {
            
            qq <- matrix(0,length(Xtrim),1,byrow=TRUE)
            
            for(i in 1:length(Xtrim)) {
              Xpoly <- (Xtrim - Xtrim[i])
              Xweight <- Kweights((Xpoly)/(BW[h]*(sd(Xtrim))*(dim(as.matrix(Xtrim))[1])^(-1/3)))
              if (colSums(Xweight!=0)<=5) {
                Xweight <- Kweights((Xpoly)/((BW[h]*(sd(Xtrim))*(dim(as.matrix(Xtrim))[1])^(-1/3))+.0001))	
              }
              qq[i] <- rq(Ytrim~cbind(Xpoly,(Xpoly^2),(Xpoly^3)), weights=Xweight, tau=tau_grid[[tau]], ci=FALSE)$coef[1.]
            }
            
            qhat[[tau]] <- qq
          #qhat[[tau]] <- 2*(Xtrim-0.5)^3 + (Xtrim-0.5)  + qnorm(tau_grid[[tau]], mean = 0, sd = 0.5)
  
            
          }
          
            
          for (del in 1:length(hpvec)) {
            
            kweight   <-  Kweights((Ptrim-deltavec[del])/hpvec[del])
            if (sum(kweight)==0)
            {
              kweight   <- Kweights((Ptrim-deltavec[del])/(hpvec[del]+0.0001))
            }
            

              
              
            CQ_stat_2 <- LocVBDNTest(Ytrim,qhat, kweight, tau_grid, B=Bsize)
            TestRes[[del]][mc,((ss-1)*3+h)] <- CQ_stat_2$teststat
            BootRes[[del]][mc,((ss-1)*3+h)] <- CQ_stat_2$teststatb

          }
          
        }   
      
    }
    
    for (del in 1:length(hpvec)) {
      
        for (h in 1:length(BW)) {
          cv1 <- quantile(BootRes[[del]][,((ss-1)*3+h)], 1-alpha1, na.rm=TRUE)
          rej1=TestRes[[del]][,((ss-1)*3+h)]>cv1
          cv2 <- quantile(BootRes[[del]][,((ss-1)*3+3)], 1-alpha2, na.rm=TRUE)
          rej2=TestRes[[del]][,((ss-1)*3+h)]>cv2
          RejRate[[del]][h,((ss-1)*2+1):((ss-1)*2+2)] <- apply(cbind(rej1,rej2),2,mean)
        }
    }
    
    
  }
  
  
  
  # Pre-allocate empty results table
  table <- lapply(1:length(hpvec), function(x) matrix(0,6,4))
  for (i in 1:length(hpvec)) {
    table[[i]] <- rbind(matrix(paste("hp/del=",hpvec[i],sep=""),1,4),matrix(paste("Gamma1=",gamma1,sep=""),1,4),matrix(paste("SS=",SS,sep=""),1,4),RejRate[[i]])
  }
  
  table = do.call(cbind,table)
  # Save csv of results
  path=paste("MC-Results-Ext-CQ",CQ,"-T1",T1,"-DGP-SE-",SE,"-spacing-Gamma1-",gamma1,"-Gamma2-",(gamma2),".csv",sep="")
  write.table(table,path, row.names=FALSE, col.names=FALSE, sep=",")
  
}

# Total Running Time:
cat(paste("Running Time: ",floor(as.matrix(proc.time()-time)[3,]/60),"mins ",round((as.matrix(proc.time()-time)[3,]-60*floor(as.matrix(proc.time()-time)[3,]/60))),"sec",sep=""))
