########################################################
## Simulation Routine: Conditional Quantile Test 1    ##
## For: Results of Table 1, Supplementary Material    ##
## Paper: Testing for Quantile Sample Selection       ##
## Authors: V. Corradi and D. Gutknecht               ##
## Date: 20-04-22                                     ##
########################################################
## Note: some of the results from Table 1 have been generated with a different seed,
## and so may numerically differ slightly. The qualitative results remain however identical.


rm(list = ls())


library(np)
library(tidyr)
library(quantreg)
library(MASS)
source("Multiplot.R")



##########################
##########################
######  FUNCTIONS   ######
##########################
##########################




#' compute Volgushev et al. (2013) test statistic
#'
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
WBootStat <- function(Bernl,n, IndX, IndZ, taueval) {
  
  
  Tn.vals  <- vector(mode="list",dim(as.matrix(taueval))[1])
  
  for (tau in 1:dim(as.matrix(taueval))[1]) {
    
      for (i in 1:dim(IndX)[2]) {
        Tn.vals[[tau]][i] <- mean((Bernl[[tau]] - taueval[tau])*IndX[,i]*IndZ[[tau]][,i])
      }
  }
  
  
  teststatb <- max(abs(sqrt(n)*unlist(Tn.vals)))
  
  return(list(teststat=teststatb))
}



#' perform Volgushev et al. (2013) test 
VBDNTest <- function(Y, IndX, IndZ, T1IndZ, qhat, taueval, B) {
  
  
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
    teststatb[b] <- WBootStat(Bernl,n , IndX, T1IndZ, taueval)$teststat
    
  }

  return(list(teststat=teststat, teststatb=teststatb))
}


##########################
##########################
### SIMULATION ROUTINE ###
##########################
##########################



time=proc.time()

# Set working directory

setwd("C:/Users/gutknech/Downloads/Simulations/Results")

# Set number of Monte Carlo replications
MC <- 999

# Set the nominal size of both tests
alpha1 <- 0.05
alpha2 <- 0.1

# Set the number of bootstrap replications (Note: Warp bootstrap - one bootstrap draw per Monte Carlo replication)
Bsize <- 1

# Choose Case 1-8 (see supplementary material, S13-S18):

case <- 1

# Define the grid of quantiles
tau_grid = c(.3,.5,.7) 

# Set the degree of selection rho={0,0.25,0.5} (Note: will be over-written for CASE VII-VIII)
rho <- 0

#Only CASE VII-VIII: Set misspecification parameter gamma={0.25,0.5} (Note: will be over-written for other cases)
gamma <- 0.25

# Set trimming constant
trim <- .025


# Set the sample sizes

SS <- matrix(c(1000, 2000),1,2)




# Results matrix
if (case<=4) {
  BW<-c(3.5,4.0,4.5)
  TestRes <- lapply(1:length(BW), function(x) matrix(0,MC,dim(SS)[2]*3))
  # Rejection Rate matrix
  RejRate <- matrix(0,length(BW),4)
  
} else if (case>4){
  hgrid<-c(1.25,1.5,1.75,2,2.25,2.5,2.75) 
  TestRes <- lapply(1, function(x) matrix(0,MC,dim(SS)[2]*3))
  # Rejection Rate matrix
  RejRate <- matrix(0,1,4)
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

## Epanechnikov (2nd order) kernel weights

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
      
      
      if (case==1) {
        Z <- rnorm(Ssize,mean=0, sd=1)
        gamma <- 0
      } else if (case==2) {
        Z <- 0.5-rbinom(Ssize, 1, 0.5)
        gamma <- 0
      } else if (case==3) {
        Z <- 1.5-rpois(Ssize,1.5)
        gamma <- 0
      } else if (case==4) {
        Z <- 3.5-ceiling(runif(Ssize, min=0, max=7))
        gamma <- 0
      }  else if (case==5) {
        Z <- rnorm(Ssize,mean=0, sd=1)
        gamma <- 0
      } else if (case==6) {
        Z <- rnorm(Ssize,mean=0, sd=1)
        gamma <- 0
      } else if (case==7) {
        Z <- rnorm(Ssize,mean=0, sd=1)
        rho <- 0
      } else if (case==8) {
        Z <- 0.5-rbinom(Ssize, 1, 0.5)
        rho <- 0
      }
      
      
      X <- runif(Ssize, min = 0, max = 1)
      E <- mvrnorm(Ssize,matrix(c(0,0),2,1),matrix(c(1,rho,rho,1),2,2,byrow=TRUE))
      S <- ((0.75*(X-0.5)+(0.75*Z))>=E[,2])
      
      
      

      Y <- (X-0.5)^3 + (X-0.5)^2 + (X-0.5) + (gamma*Z)  + 0.5*E[,1]
        

      
      # Construct the propensity score
      if (case==1) {
        Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)
      } else if (case==2) {
        Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)
      } else if (case==3) { 
        Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)      
      } else if (case==4) {
        Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)
      } else if (case==5) {
        Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)
      } else if (case==6) {
        bw = npregbw(formula = as.numeric(S) ~ factor(Z) + X,nmulti=1)$bw
        fZhat <- c(np::npksum(txdat=cbind(Z,X), leave.one.out=FALSE,  bws=bw, ckertype="epanechnikov")$ksum)     
        Pz <- c(np::npksum(txdat=cbind(Z,X), tydat=as.numeric(S), leave.one.out=FALSE, bws=bw, ckertype="epanechnikov")$ksum) / (fZhat+0.00001) 
        rm(fZhat)
      }  else if (case==7) {
        Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)
      }  else if (case==8) {
        Pz <- pnorm((0.75*(X-0.5)+(0.75*Z)), 0, 1)
      }    
    
      
      Y <- Y[S]
      X <- X[S]
      Pz <- Pz[S]

      filter <- (X<quantile(X,probs=trim))|(X>quantile(X,probs=1-trim))
      Xtrim <- X[!filter]
      Ptrim <- Pz[!filter]
      Ytrim <- Y[!filter]
      rm(filter)
      
      
      #########################
      #  Computation of Test  #
      #########################
      if (case<=4) {  
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

            
          }
          
          
          

            epsilonhat <- vector(mode="list",dim(as.matrix(tau_grid))[1])
            CFp <- vector(mode="list",dim(as.matrix(tau_grid))[1])
            T1IndZ <- vector(mode="list",dim(as.matrix(tau_grid))[1])
            
            
            for (tau in 1:dim(as.matrix(tau_grid))[1]) {
              epsilonhat[[tau]] <- (Ytrim - qhat[[tau]])
              Ptrim <- as.matrix(Ptrim)
              bwFhat <- (2.2*(apply(cbind(Xtrim,as.numeric(epsilonhat[[tau]])),2,sd))*(dim(as.matrix(Ptrim))[1])^(-1/6))
             fphat <- np::npksum(txdat=data.frame(Xtrim,epsilonhat[[tau]]),exdat=data.frame(Xtrim,matrix(0,length(Ptrim),1)), leave.one.out=FALSE, bws=bwFhat, ckertype="epanechnikov")$ksum
              
              for (t in 1:length(Ptrim)) {
                CFp[[tau]][t] <-  c(np::npksum(tydat=as.numeric(Ptrim<= Ptrim[t]),txdat=data.frame(Xtrim,epsilonhat[[tau]]),exdat=data.frame(Xtrim[t],0), leave.one.out=FALSE,  bws=bwFhat, ckertype="epanechnikov")$ksum)/fphat[t]
              }
              T1IndZ[[tau]] <- T1IndMat(Ptrim, CFp[[tau]])
            }

          
          IndX <- IndMat(Xtrim)
          IndZ <- IndMat(Ptrim)
          
          CQ_stat_1 = VBDNTest(Ytrim,IndX,IndZ,T1IndZ,qhat, tau_grid, B=Bsize)

          
          TestRes[[h]][mc,((ss-1)*3+2)] <- CQ_stat_1$teststat
          TestRes[[h]][mc,((ss-1)*3+3)] <- CQ_stat_1$teststatb
        }
      }  else if (case>4) {
        
        qhat <- vector(mode="list",dim(as.matrix(tau_grid))[1])
        n <- dim(as.matrix(Xtrim))[1]
        for (tau in 1:dim(as.matrix(tau_grid))[1]) {
          
          
          n <- dim(as.matrix(Xtrim))[1]
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
        
        

          
          epsilonhat <- vector(mode="list",dim(as.matrix(tau_grid))[1])
          CFp <- vector(mode="list",dim(as.matrix(tau_grid))[1])
          T1IndZ <- vector(mode="list",dim(as.matrix(tau_grid))[1])
          
          
          for (tau in 1:dim(as.matrix(tau_grid))[1]) {
            epsilonhat[[tau]] <- (Ytrim - qhat[[tau]])
            Ptrim <- as.matrix(Ptrim)
            hCV = mnpCV(hgrid,Ptrim,cbind(Xtrim,as.numeric(epsilonhat[[tau]])))
            bwFhat <- (hCV*(apply(cbind(Xtrim,as.numeric(epsilonhat[[tau]])),2,sd))*(dim(as.matrix(Ptrim))[1])^(-1/6))
            #bwFhat <- (2.2*(apply(cbind(Xtrim,as.numeric(epsilonhat[[tau]])),2,sd))*(dim(as.matrix(Ptrim))[1])^(-1/6))
            fphat <- np::npksum(txdat=data.frame(Xtrim,epsilonhat[[tau]]),exdat=data.frame(Xtrim,matrix(0,length(Ptrim),1)), leave.one.out=FALSE, bandwidth.divide=TRUE, bws=bwFhat, ckertype="epanechnikov")$ksum
            
            for (t in 1:length(Ptrim)) {
              CFp[[tau]][t] <-  c(np::npksum(tydat=as.numeric(Ptrim<= Ptrim[t]),txdat=data.frame(Xtrim,epsilonhat[[tau]]),exdat=data.frame(Xtrim[t],0), leave.one.out=FALSE, bandwidth.divide=TRUE, bws=bwFhat, ckertype="epanechnikov")$ksum)/fphat[t]
            }
            T1IndZ[[tau]] <- T1IndMat(Ptrim, CFp[[tau]])
          }
        
        
        
        IndX <- IndMat(Xtrim)
        IndZ <- IndMat(Ptrim)
        
        CQ_stat_1 = VBDNTest(Ytrim,IndX,IndZ,T1IndZ,qhat, tau_grid, B=Bsize)

        
        TestRes[[1]][mc,((ss-1)*3+2)] <- CQ_stat_1$teststat
        TestRes[[1]][mc,((ss-1)*3+3)] <- CQ_stat_1$teststatb
      }    
      
      
    }
    
    if (case<=4) {  
      
      for (h in 1:length(BW)) {
        
        cv1 <- quantile(TestRes[[h]][,((ss-1)*3+3)], 1-alpha1, na.rm=TRUE)
        rej1=TestRes[[h]][,((ss-1)*3+2)]>cv1
        cv2 <- quantile(TestRes[[h]][,((ss-1)*3+3)], 1-alpha2, na.rm=TRUE)
        rej2=TestRes[[h]][,((ss-1)*3+2)]>cv2
        RejRate[h,((ss-1)*2+1):((ss-1)*2+2)] <- apply(cbind(rej1,rej2),2,mean)
        
      }
    } else if (case>4) {
      cv1 <- quantile(TestRes[[1]][,((ss-1)*3+3)], 1-alpha1, na.rm=TRUE)
      rej1=TestRes[[1]][,((ss-1)*3+2)]>cv1
      cv2 <- quantile(TestRes[[1]][,((ss-1)*3+3)], 1-alpha2, na.rm=TRUE)
      rej2=TestRes[[1]][,((ss-1)*3+2)]>cv2
      RejRate[1,((ss-1)*2+1):((ss-1)*2+2)] <- apply(cbind(rej1,rej2),2,mean)
      
    }
    
    
    print(RejRate[,ss])
  }
  
  
  # Pre-allocate empty results table
  table=matrix(NA,0,ncol(RejRate))
  
  table=rbind(matrix(paste("rho=",rho,sep=""),1,4),matrix(paste("SS=",kronecker(SS,matrix(1,1,2)),sep=""),1,4),RejRate)
  
  # Save csv of results
  path=paste("MC-Results-CQ-T1-DGP-CASE-",case,"-Rho-",(rho*4),".csv",sep="")
  write.table(table,path, row.names=FALSE, col.names=FALSE, sep=",")
 

# Total Running Time:
cat(paste("Running Time: ",floor(as.matrix(proc.time()-time)[3,]/60),"mins ",round((as.matrix(proc.time()-time)[3,]-60*floor(as.matrix(proc.time()-time)[3,]/60))),"sec",sep=""))
