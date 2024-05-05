########################################################
## Simulation Routine: Conditional Quantile Test 2    ##
## For: Results Tables 2 & 3, Supplementary Material  ##
## Paper: Testing for Quantile Sample Selection       ##
## Authors: V. Corradi and D. Gutknecht               ##
## Date: 08-09-22                                     ##
########################################################

rm(list = ls())

## Note: working directory has to be changed to the user's working directory. ALL OUTPUT WILL BE SAVED IN WORKING DIRECTORY
setwd("C:/Users/gutknech/Downloads/Simulations/Results")


library(np)
library(tidyr)
library(quantreg)
library(MASS)


##########################
##########################
######  FUNCTIONS   ######
##########################
##########################

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

# Bootstrap statistic

LocWBootStat <- function(Bernl,kweights , taueval) {
  
  Tn.vals <- vector(mode="list",dim(as.matrix(taueval))[1])
  Var.vals <- vector(mode="list",dim(as.matrix(taueval))[1])
  
  for (tau in 1:dim(as.matrix(taueval))[1]) {
    
    Tn.vals[[tau]] <- sum((Bernl[[tau]] - taueval[tau])*kweights)
    Var.vals[[tau]] <- mean((Bernl[[tau]] - taueval[tau])^2)*sum(kweights)*(3/5)          
  } 
  # bootstrap statistic
  
  teststatb <- max(abs(unlist(Tn.vals)/sqrt(unlist(Var.vals))),na.rm=TRUE)
  
  return(list(teststat=teststatb))
}


# Final Function

LocVBDNTest <- function(Y, qhat, kweights, taueval, B) {
    
  
  n <- length(Y)

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

## Kernel Weights

Kweights <- function(arg) {
  
  arg <- as.matrix(arg)
  return((3/4)*(1-(arg^2))*(as.numeric(abs(arg)<=1)))
  
}



##########################
##########################
### SIMULATION ROUTINE ###
##########################
##########################


time=proc.time()

# Set number of Monte Carlo replications
MC <- 999

# Set the nominal size of both tests
alpha1 <- 0.05
alpha2 <- 0.1


# Set the number of bootstrap replications (Note: Warp bootstrap - one bootstrap draw per Monte Carlo replication)
Bsize <- 1



# Set the degree of selection (correlation parameter fixed at 0.25 in Tables 2 and 3, Supplementary Material)    
rho <- 0.25

  

# Define the grid of quantiles
tau_grid = c(.3,.5,.7) 


# Set trimming constant
trim <- .025



# Set the sample sizes

SS <- matrix(c(1000, 2000),1,2)


##########################
##########################
### SIMULATION ROUTINE ###
##########################
##########################

# Loop over Case I*-VI* (see supplementary material, S13-S18):
# CASE I* - 1; CASE II* - 2; CASE III* - 3; CASE IV* - 4; CASE V* - 5; CASE VI* - 6; 
  
for (casestar in 1:6) {  
  
  for (BWsel in 1:2) {  
  # Set random number seed and get start time
  
    
  if (BWsel==1) {
    
    # Scaling Factors (cf. Section S2, supplementary material)
    BW<-c(3.5,4.0,4.5)
    
    if (casestar<=3) {
      
      deltavec <- c(0.95,0.95,0.975,.98)
      hpvec    <- c(0.075,0.05,0.03,.025)
      gamma2 <- 1
      
    } else if (casestar>3){
      
      deltavec <- c(0.95,0.95)
      hpvec    <- c(0.075,0.05)
      gamma2 <- 0.5
      
    }
  
    # Rejection Rate list
    
    RejRate <- lapply(1:length(hpvec), function(x) matrix(0,length(BW),4))
    
      
    TestRes <- lapply(1:length(hpvec), function(x) matrix(0,MC,dim(SS)[2]*3))
    BootRes <- lapply(1:length(hpvec), function(x) matrix(0,MC,dim(SS)[2]*3))
    
    # set random seed number
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
      
        
        if (casestar<=3) {
          Z <- rnorm(Ssize,mean=0, sd=1)
        } else if (casestar>3) {
          Z <- 0.5-rbinom(Ssize, 1, 0.5)
        } 
        
      
      if (casestar==1|casestar==4) {
        gamma1 <- 0
      }
      
      if (casestar==2|casestar==5) {
        gamma1 <- 0.25
      }
      
      if (casestar== 3|casestar==6) {
        gamma1 <- 0.5
      }
      
      
        X <- runif(Ssize, min = 0, max = 1)
        E <- mvrnorm(Ssize,matrix(c(0,0),2,1),matrix(c(1,rho,rho,1),2,2,byrow=TRUE))
        S <- ((0.75*(X-0.5)+(0.75*Z))>=(gamma2*E[,2]))
        
        
        Y <- (X-0.5)^3 + (X-0.5)^2 + (X-0.5)+ (gamma1*Z)  + 0.5*E[,1]
          
        Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, (gamma2*1))

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
                Xweight <- Kweights((Xpoly)/((BW[h]*(sd(Xtrim))*(dim(as.matrix(Xtrim))[1])^(-1/3))+.1))	
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
    table[[i]] <- rbind(matrix(paste("hp/del=",hpvec[i],sep=""),1,4),matrix(paste("Gamma1=",gamma1,sep=""),1,4),matrix(paste("SS=",kronecker(SS,matrix(1,2)),sep=""),1,4),RejRate[[i]])
  }
  
  table = do.call(cbind,table)
  # Save csv of results
  path=paste("MC-Results-NO-BW-DGP-CASESTAR-",casestar,".csv",sep="")
  write.table(table,path, row.names=FALSE, col.names=FALSE, sep=",")

  ##########################################################
  # File saved as:                                         #
  # MC-Results-NO-BW-DGP-CASESTAR-XXX.csv                  #
  # XXX: CASE I*-VI* in roman letters                      #
  ##########################################################
  # Format of the saved csv result file (CASE I*-III*):    #
  # Rows:                                                  #
  # - Row 1: hp/ del value                                 #
  # - Row 2: value of Gamma1  (misspecification)           #
  # - Row 3: sample size n ("SS")                          #
  # - Row 4: results for c=3.5                             #  
  # - Row 5: results for c=4.0                             #
  # - Row 6: results for c=4.5                             #
  # Columns:                                               #                                       #
  # - Column 1: n=1000 and alpha=0.05 and del=0.95/hp=0.75 #
  # - Column 2: n=1000 and alpha=0.10 and del=0.95/hp=0.75 #
  # - Column 3: n=2000 and alpha=0.05 and del=0.95/hp=0.75 #
  # - Column 4: n=2000 and alpha=0.10 and del=0.95/hp=0.75 #
  # - Column 5: n=1000 and alpha=0.05 and del=0.95/hp=0.05 #
  # - Column 6: n=1000 and alpha=0.10 and del=0.95/hp=0.05 #
  # - Column 7: n=2000 and alpha=0.05 and del=0.95/hp=0.05 #
  # - Column 8: n=2000 and alpha=0.10 and del=0.95/hp=0.05 #
  # - Column 9: n=1000 and alpha=0.05 and del=0.975/hp=0.03#
  # - Column10: n=1000 and alpha=0.10 and del=0.975/hp=0.03#
  # - Column11: n=2000 and alpha=0.05 and del=0.975/hp=0.03#
  # - Column12: n=2000 and alpha=0.10 and del=0.975/hp=0.03#
  # - Column13: n=1000 and alpha=0.05 and del=0.98/hp=0.02 #
  # - Column14: n=1000 and alpha=0.10 and del=0.98/hp=0.02 #
  # - Column15: n=2000 and alpha=0.05 and del=0.98/hp=0.02 #
  # - Column16: n=2000 and alpha=0.10 and del=0.98/hp=0.02 #
  ##########################################################
  # Format of the saved csv result file (CASE IV*-VI*):    #
  # Rows:                                                  #
  # - Row 1: hp/ del value                                 #
  # - Row 2: value of Gamma1 (misspecification)            #
  # - Row 3: sample size n ("SS")                          #
  # - Row 4: results for c=3.5                             #  
  # - Row 5: results for c=4.0                             #
  # - Row 6: results for c=4.5                             #
  # Columns:                                               #                                       #
  # - Column 1: n=1000 and alpha=0.05 and del=0.95/hp=0.75 #
  # - Column 2: n=1000 and alpha=0.10 and del=0.95/hp=0.75 #
  # - Column 3: n=2000 and alpha=0.05 and del=0.95/hp=0.75 #
  # - Column 4: n=2000 and alpha=0.10 and del=0.95/hp=0.75 #
  # - Column 5: n=1000 and alpha=0.05 and del=0.95/hp=0.05 #
  # - Column 6: n=1000 and alpha=0.10 and del=0.95/hp=0.05 #
  # - Column 7: n=2000 and alpha=0.05 and del=0.95/hp=0.05 #
  # - Column 8: n=2000 and alpha=0.10 and del=0.95/hp=0.05 #
  ##########################################################
  
  }
    
  if (BWsel==2) {
    
    
    if (casestar<=3) {
      
      gamma2 <- 1
      
    } else if (casestar>3){
      
      gamma2 <- 0.5
      
    }
    
    # Rejection Rate list
    RejRate <- matrix(0,2,2)
    
    
    
    TestRes <- matrix(0,MC,dim(SS)[2])
    BootRes <- matrix(0,MC,dim(SS)[2])
    
    hgrid<-c(1.25,1.5,1.75,2,2.25,2.5,2.75) 
    
    
    etavec <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
    eps <- 0.1
    
    # Set random number seed
    
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
        
        
        
        
        if (casestar<=3) {
          Z <- rnorm(Ssize,mean=0, sd=1)
        } else if (casestar>3) {
          Z <- 0.5-rbinom(Ssize, 1, 0.5)
        } 
        
        if (casestar==1|casestar==4) {
          gamma1 <- 0
        }
        
        if (casestar==2|casestar==5) {
          gamma1 <- 0.25
        }
        
        if (casestar== 3|casestar==6) {
          gamma1 <- 0.5
        }
        
        X <- runif(Ssize, min = 0, max = 1)
        E <- mvrnorm(Ssize,matrix(c(0,0),2,1),matrix(c(1,rho,rho,1),2,2,byrow=TRUE))
        S <- ((0.75*(X-0.5)+(0.75*Z))>=(gamma2*E[,2]))
        
        
        Y <- (X-0.5)^3 + (X-0.5)^2 + (X-0.5)+ (gamma1*Z)  + 0.5*E[,1]
        
        Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, (gamma2*1))
        
        Y <- Y[S]
        X <- X[S]
        Pz <- Pz[S]
        
        filter <- (X<quantile(X,probs=trim))|(X>quantile(X,probs=1-trim))
        Xtrim <- X[!filter]
        Ptrim <- Pz[!filter]
        Ytrim <- Y[!filter]
        rm(filter)
        
        
        ########################
        # Conditional Quantile #
        ########################
        
        qhat <- vector(mode="list",dim(as.matrix(tau_grid))[1])
        n <- dim(as.matrix(Xtrim))[1]
        for (tau in 1:dim(as.matrix(tau_grid))[1]) {
          
          
          CV_err_h =rep(NA,length(hgrid))
          for(j in 1:length(hgrid)){
            h_using = (hgrid[j]*sd(Xtrim)*(n^(-1/5)))
            CV_err =rep(NA, n)
            qq <- matrix(0,(n-1),1,byrow=TRUE)
            for(i in 1:n){
              X_val = Xtrim[i]
              Y_val = Ytrim[i]
              X_tr = Xtrim[-i]
              Y_tr = Ytrim[-i]
              Xweight <- Kweights((X_tr-X_val)/h_using)
              if (colSums(Xweight!=0)<=5) {
                Xweight <- Kweights((X_tr-X_val)/(h_using+.05))
              }
              Xpoly <- (X_tr-X_val)
              qq[i] <- rq(Y_tr~Xpoly, weights=Xweight, tau=tau_grid[[tau]], ci=FALSE)$coef[1.]
              CV_err[i] = (Y_val - qq[i])*(tau_grid[[tau]]-((Y_val - qq[i])<=0))
              
            }
            CV_err_h[j] =mean(CV_err)
          }
          
          BW <- (hgrid[which(CV_err_h ==min(CV_err_h))]*sd(Xtrim)*(n^(-1/5)))
          
          qq <- matrix(0,length(Xtrim),1,byrow=TRUE)
          
          for(i in 1:length(Xtrim)) {
            
            Xweight <- Kweights((Xtrim-Xtrim[i])/BW)
            if (colSums(Xweight!=0)<=5) {
              Xweight <- Kweights((Xtrim-Xtrim[i])/(BW+.05))
            }
            Xpoly <- (Xtrim - Xtrim[i])
            qq[i] <- rq(Ytrim~cbind(Xpoly,(Xpoly^2),(Xpoly^3)), weights=Xweight, tau=tau_grid[[tau]], ci=FALSE)$coef[1.]
          }
          
          qhat[[tau]] <- qq
          
        }
        kweightfinal <- matrix(NA,1,length(Ptrim))      
        etaind <- length(etavec)
        while (etaind > 0) {
          eta <- etavec[etaind]      
          hp <- log(n)*(n^(-((1+eps)/(1+eta+eps))))
          delta <- 1-(hp^(1/(1+eps)))
          kweight <- Kweights((Ptrim-delta)/hp)
          if ((sum(kweight)/(n*hp))>=0.1) {
            kweightfinal <- Kweights((Ptrim-delta)/hp)
          } else { break }
          etaind = etaind-1
        }
        
        CQ_stat_2 <- LocVBDNTest(Ytrim,qhat, kweightfinal, tau_grid, B=Bsize)
        TestRes[mc,ss] <- CQ_stat_2$teststat
        BootRes[mc,ss] <- CQ_stat_2$teststatb
        
        
        
        
        
      }
      
      cv1 <- quantile(BootRes[,ss], 1-alpha1, na.rm=TRUE)
      rej1=TestRes[,ss]>cv1
      cv2 <- quantile(BootRes[,ss], 1-alpha2, na.rm=TRUE)
      rej2=TestRes[,ss]>cv2
      RejRate[,ss] <- apply(cbind(rej1,rej2),2,mean)
      
      
    }
    
    
    
    # Pre-allocate empty results table
    table <- matrix(0,1,2)
    table <- rbind(matrix(paste("hp-hat",hp,sep=""),1,2),matrix(paste("Gamma1=",gamma1,sep=""),1,2),matrix(paste("SS=",SS,sep=""),1,2),RejRate)
    # Save csv of results
    path=paste("MC-Results-BW-Select--DGP-CASESTAR-",casestar,".csv",sep="")
    write.table(table,path, row.names=FALSE, col.names=FALSE, sep=",")
    
    #########################################################
    # File saved as:                                        #
    # MC-Results-BW-Select--DGP-CASESTAR-XXX.csv            #
    #                                                       #
    # Format of the saved csv result file:                  #
    # Rows:                                                 #
    # - Row 1: data driven bandwidth value of last run      #
    # - Row 2: value of Gamma1 (misspecification)           #
    # - Row 3: sample size n ("SS")                         #
    # - Row 4: results for alpha=0.05                       #  
    # - Row 5: results for alpha=0.10                       #
    # Columns:                                              #
    # - Column 1: n=1000                                    #
    # - Column 2: n=2000                                    #
    #########################################################
    
  }
    
    
  }
  
}

# Total Running Time:
cat(paste("Running Time: ",floor(as.matrix(proc.time()-time)[3,]/60),"mins ",round((as.matrix(proc.time()-time)[3,]-60*floor(as.matrix(proc.time()-time)[3,]/60))),"sec",sep=""))
