###########################################
## Simulation Routine - Conditional Mean ##
###########################################

## Note: this routine requires executing "CMTest1FINAL.R" and "CMTest2FINAL.R" first.

rm(list = ls())

# load packages
library(np)
library(tidyr)
library(quantreg)
library(MASS)
source("Multiplot.R")


time=proc.time()

# Set working directory

setwd("C:/Users/gutknech/Downloads/Simulations/Results")

# Set number of Monte Carlo replications
MC <- 999

# Set the nominal size of both tests
alpha1 <- 0.05
alpha2 <- 0.1

# Set the number of bootstrap replications (Note: "Warp speed" bootstrap is used, see Giacomini, Politis, and White (2013))
Bsize <- 1

# Choose the global or local test (T1==1: `global'  (First Test; Section 3, pp.7-12); T1==0: `local' (Second Test; Section 4, pp.12-18)
T1 <- 1

# Set parameters of Second (Local) Test:
if (T1==0) {
  
  deltavec <- c(0.95,0.975,0.98)
  hpvec    <- c(0.05,0.025,0.02)
  
  
}

# Define the grid of quantile levels
tau_grid = c(.3,.4,.5,.6,.7) 


# Set the degree of selection 
rho <- 0.25

# Set trimming constant (trims away the outer "trim"% from the observations on X)
trim <- .025


# Cond. Mean Bandwidths:
BW<-c(0.125,0.25,0.5)



# Set the sample sizes
SS <- matrix(c(400, 1000),1,2)




# Rejection Rate matrix
RejRate <- matrix(0,length(BW),4)

# Results matrix
TestRes <- lapply(1:length(BW), function(x) matrix(0,MC,dim(SS)[2]*3))


## Evaluation matrix (indicator) functions:

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


# Second Order Epanechnikov  kernel weight function (required for second test only)
Kweights <- function(arg) {
  
  arg <- as.matrix(arg)
  return((3/4)*(1-(arg^2))*(as.numeric(abs(arg)<=1)))
  
}









##########################
##########################
### SIMULATION ROUTINE ###
##########################
##########################

for (SE in 1:6) {
# Set random number seed
  
  set.seed(12345)
  
for (ss in 1:dim(SS)[2]) {
  
  
  Ssize <- SS[ss]

  # Set bandwidth for second (local) test:
  hp <- Ssize^(-(1/5))

  
  for (mc in 1:MC) {
    ########################
    # Data Generation      #
    ########################
    
    # Output simulation rounds:
    if ((mc %in% seq(100,100*floor(MC/100),100))==TRUE){
      cat(mc)
    } else {
      cat('.')
    }
    # Select nature of instrument (continuous, poisson, binomial):
    if (SE==1) {
      Z <- rnorm(Ssize,mean=0, sd=1)
    } else if (SE==2) {
      Z <- 0.5-rbinom(Ssize, 1, 0.5)
    }  else if (SE==3) {
      Z <- 1.5-rpois(Ssize,1.5)
    } else if (SE==4) {
      Z <- 3.5-ceiling(runif(Ssize, min=0, max=7))
    } else if (SE==5) {
      Z <- rnorm(Ssize,mean=0, sd=1)
    } else if (SE==6) {
      Z <- 0.5-rbinom(Ssize, 1, 0.5)
    } 
    
    X <- runif(Ssize, min = 0, max = 1)
    E <- mvrnorm(Ssize,matrix(c(0,0),2,1),matrix(c(1,rho,rho,1),2,2,byrow=TRUE))
    S <- ((0.75*(X-0.5)+(0.75*Z))>=E[,2])
    
    
    

      Y <- X^2 + 0.5*X + 0.5*E[,1]
      
    # Construct the propensity score
    if (SE==1) {
      Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)
    } else if (SE==2) {
      Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)
    } else if (SE==3) { 
      
      Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)      
      
    } else if (SE==4) {
      Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)
      
    }  else if (SE==5) {
      bw = npregbw(formula = as.numeric(S) ~ Z + X,nmulti=1)$bw
      fZhat <- c(np::npksum(txdat=cbind(Z,X), leave.one.out=FALSE,  bws=bw, ckertype="epanechnikov")$ksum)     
      Pz <- c(np::npksum(txdat=cbind(Z,X), tydat=as.numeric(S), leave.one.out=FALSE, bws=bw, ckertype="epanechnikov")$ksum) / (fZhat+0.00001) 
      rm(fZhat)
    }  else if (SE==6) {
      bw = npregbw(formula = as.numeric(S) ~ factor(Z) + X,nmulti=1)$bw
      fZhat <- c(np::npksum(txdat=cbind(Z,X), leave.one.out=FALSE,  bws=bw, ckertype="epanechnikov")$ksum)     
      Pz <- c(np::npksum(txdat=cbind(Z,X), tydat=as.numeric(S), leave.one.out=FALSE, bws=bw, ckertype="epanechnikov")$ksum) / (fZhat+0.00001) 
      rm(fZhat)
    }    
    
    Y <- Y[S]
    X <- X[S]
    Pz <- Pz[S]
    
    # Generate selected sample:  
    
    filter <- (X<quantile(X,probs=trim))|(X>quantile(X,probs=1-trim))
    Xtrim <- X[!filter]
    Ptrim <- Pz[!filter]
    Ytrim <- Y[!filter]
    rm(filter)
    

      
      #########
      # Tests #
      #########
      
        
        if (T1==1) {           
          ## FIRST TEST

          for (h in 1:length(BW)) {
            
          
          IndX <- IndMat(Xtrim)
          IndZ <- IndMat(Ptrim)
          
          CM_stat_1 = DMTest(Ytrim,Xtrim,IndX,IndZ, size=alpha1, B=Bsize, a=(BW[h]*(dim(as.matrix(Ytrim))[1])^(-1/3)), ckertype="epanechnikov",stat="KS") 
          TestRes[[h]][mc,((ss-1)*3+2)] <- CM_stat_1$teststat
          TestRes[[h]][mc,((ss-1)*3+3)] <- CM_stat_1$teststatb
          
          }  
          
        } else if (T1==0) {
          
          
          kweight   <-  Kweights((Ptrim-delta)/hp)
          if (sum(kweight)==0)
          {
            kweight   <- Kweights((Ptrim-delta)/(hp+0.01))
            if (sum(kweight)==0) 
            {
              kweight   <- Kweights((Ptrim-delta)/(abs(max(Ptrim)-delta)+0.01))  
            }
          }
          
          n <- length(kweight)
          
          for (h in 1:length(BW)) {
            
            
          # nonparametric regression of Y on X
          fhat <- c(np::npksum(txdat=Xtrim, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=BW[h], ckertype="epanechnikov")$ksum) / (n*BW[h])
          Yhat <- c(np::npksum(txdat=Xtrim, tydat=Ytrim, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=BW[h], ckertype="epanechnikov")$ksum) / (n*BW[h]*fhat)
          epsilonhat <- Ytrim-Yhat
          Testvar <- epsilonhat*fhat
              
          ## SECOND TEST
          
          
#          CM_stat_2 = LocDMTest(Testvar,kweight, size=alpha1, B=Bsize)
          CM_stat_2 = LocDMTest(Testvar,kweight, size=alpha1) 
          
          TestRes[[h]][mc,((ss-1)*3+1)] <- as.numeric(CM_stat_2$rej)
          TestRes[[h]][mc,((ss-1)*3+2)] <- as.numeric(CM_stat_2$teststat)
          TestRes[[h]][mc,((ss-1)*3+3)] <- as.numeric(CM_stat_2$cv)
          
          }
          
        }
        
      
      
     
    
  }
  
  for (h in 1:length(BW)) {
    
    cv1 <- quantile(TestRes[[h]][,((ss-1)*3+3)], 1-alpha1, na.rm=TRUE)
    rej1=TestRes[[h]][,((ss-1)*3+2)]>cv1
    cv2 <- quantile(TestRes[[h]][,((ss-1)*3+3)], 1-alpha2, na.rm=TRUE)
    rej2=TestRes[[h]][,((ss-1)*3+2)]>cv2
    RejRate[h,((ss-1)*2+1):((ss-1)*2+2)] <- apply(cbind(rej1,rej2),2,mean)
    
  }


	print(RejRate[,ss])
}

  
  # Pre-allocate empty results table
  table=matrix(NA,0,ncol(RejRate))
  
  table=rbind(matrix(paste("rho=",rho,sep=""),1,4),matrix(paste("SS=",SS,sep=""),1,4),RejRate)
  
  # Save csv of results
  path=paste("MC-Results-CM-DGP-SE-",SE,"-Rho-",(rho*4),".csv",sep="")
  write.table(table,path, row.names=FALSE, col.names=FALSE, sep=",")
  
}

# Total Running Time:
cat(paste("Running Time: ",floor(as.matrix(proc.time()-time)[3,]/60),"mins ",round((as.matrix(proc.time()-time)[3,]-60*floor(as.matrix(proc.time()-time)[3,]/60))),"sec",sep=""))
