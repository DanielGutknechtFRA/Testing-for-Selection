########################################################
## Application Routine                                ##
## For: Results of Table 1 and 2, 1995-1997 Subsample ##
## Paper: Testing for Quantile Sample Selection       ##
## Authors: V. Corradi and D. Gutknecht               ##
## Date: 08-09-22                                     ##
########################################################

rm(list = ls())

## Note: working directory has to be changed to the user's working directory. ALL OUTPUT WILL BE SAVED IN WORKING DIRECTORY
setwd("C:/Users/gutknecht/Mannheim/PROJECTS/CorradiGutknecht/Data/AB-Data/14030_Data_and_Programs/Data_files")



###################
##   FUNCTIONS   ##
##  First Test   ##
###################


#' compute (global) quantile test statistic (cf. Section 3)
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
VBDNTest <- function(Y, IndX, IndZ, T1IndZ, qhat, taueval, size, B) {
  
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
    teststatb[b] <- WBootStat(Bernl,n , IndX, T1IndZ, taueval)$teststat
    
  }
  cv <- quantile(teststatb, 1-size, na.rm=TRUE)
  
  Fn <- ecdf(teststatb)
  pval <- 1-Fn(teststat)
  
  return(list(teststat=teststat, cv=cv, rej=teststat>cv, pval=pval,n=length(Y)))
}

###################
##   FUNCTIONS   ##
##  Second Test  ##
###################


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
  
  return(list(teststat=teststat, cv=cv, rej=teststat>cv, pval=pval,n=sum(kweights>0)))
}


## Evaluation matrix functions
IndMatX <- function(X) {
  X <- as.matrix(X)
  indx <- matrix(NA,dim(X)[1],dim(X)[1])
  for( i in 1:dim(X)[1]) {
    indx[,i] <- apply((X<=X[i,]), 1, prod) 
  }
  return(indx)
}

IndMatZ <- function(Z) {
  Z <- as.matrix(Z)
  indz <- matrix(NA,dim(Z)[1],dim(Z)[1])
  for (i in 1:dim(Z)[1]) {
    indz[,i] <- as.numeric(Z<=Z[i])
  }
  return(indz)
}

T1IndMatZ <- function(Z,Fz) {
  Z <- as.matrix(Z)
  indz <- matrix(NA,dim(Z)[1],dim(Z)[1])
  for (i in 1:dim(Z)[1]) {
    indz[,i] <- ((Z<=Z[i])-Fz[i])
  }
  return(indz)
}



################
################
## Begin File ##
################
################


library(np)
library(tidyr)
library(quantreg)
# set seed
set.seed(12345)

time=proc.time()


## Define tau grid:
tau_grid = c(.1,.2,.3,.4,.5,.6,.7,.8,.9)  



# Males
data1=read.csv("data_1.csv")
# Females
data2=read.csv("data_2.csv")
# Read in bandwidth used in the paper:
BW <- read.csv(file = "BWs-1995-1997.csv",header =TRUE)
bwPf <- BW$bwPf
bwPm <- BW$bwPm
bwf <- BW$bwf[!is.na(BW$bwf)]
bwm <- BW$bwm[!is.na(BW$bwm)]
bwCFf <- BW$bwCFf[!is.na(BW$bwCFf)]
bwCFm <- BW$bwCFm[!is.na(BW$bwCFm)]




# Males == 0
data1 = cbind(matrix(0,dim(data1)[1],1),data1)
colnames(data1)[1] = c("gender")
# Females == 1
data2 = cbind(matrix(1,dim(data2)[1],1),data2)
colnames(data2)[1] = c("gender")
data = rbind(data1,data2)

datalist = list()

# Keep observations between 1995 and 1997
yrfilter =  ((data$year >= 95) & (data$year<=97))
datatmp  =  data[yrfilter,]
rm(yrfilter)
datatmp = cbind(datatmp,matrix(0,dim(datatmp)[1],3))

filter95 = (datatmp$year == 95)
filter96 = (datatmp$year == 96)
filter97 = (datatmp$year == 97)
datatmp[filter95,(dim(datatmp)[2]-2)] = 1
datatmp[filter96,(dim(datatmp)[2]-1)] = 1
datatmp[filter97,(dim(datatmp)[2])] = 1
rm(filter95,filter96,filter97)


dataFL = datatmp

colnames(dataFL)[(dim(dataFL)[2]-2):dim(dataFL)[2]] <- c("yr_95", "yr_96", "yr_97")
dataFL = data.frame(lw=dataFL$lw, work = dataFL$work,  benefit=dataFL$benefit , gender=factor(dataFL$gender), married=factor(dataFL$married), ed17=factor(dataFL$ed17),ed18=factor(dataFL$ed18), age=ordered(dataFL$age), reg_d1=factor(dataFL$reg_d1), reg_d2=factor(dataFL$reg_d2), reg_d3=factor(dataFL$reg_d3),reg_d4=factor(dataFL$reg_d4),reg_d5=factor(dataFL$reg_d5),reg_d6=factor(dataFL$reg_d6),reg_d7=factor(dataFL$reg_d7),reg_d8=factor(dataFL$reg_d8),reg_d9=factor(dataFL$reg_d9),reg_d10=factor(dataFL$reg_d10),reg_d11=factor(dataFL$reg_d11),kids_d1=factor(dataFL$kids_d1),kids_d2=factor(dataFL$kids_d2),kids_d3=factor(dataFL$kids_d3),kids_d4=factor(dataFL$kids_d4),kids_d5=factor(dataFL$kids_d5),kids_d6= factor(dataFL$kids_d6),yr_95= factor(dataFL$yr_95),yr_96= factor(dataFL$yr_96),yr_97= factor(dataFL$yr_97))
genderfil <- (dataFL$gender==1)
dataFLf <- dataFL[genderfil,]
dataFLm <- dataFL[!genderfil,]
rm(genderfil)

##################################################
## Construct Pieces Required for Test Statistic ##
##################################################



  ## Construct Nonparametric Propensity score
  
  pZf =  predict(npreg(formula = as.numeric(dataFLf$work) ~ dataFLf$benefit + dataFLf$married  + dataFLf$ed17 + dataFLf$ed18 + dataFLf$age  + dataFLf$reg_d1 + dataFLf$reg_d2 + dataFLf$reg_d3 + dataFLf$reg_d4 + dataFLf$reg_d5 + dataFLf$reg_d6 + dataFLf$reg_d7 + dataFLf$reg_d8 + dataFLf$reg_d9 + dataFLf$reg_d10 + dataFLf$reg_d11  + dataFLf$kids_d1 + dataFLf$kids_d2 + dataFLf$kids_d3 + dataFLf$kids_d4 + dataFLf$kids_d5 + dataFLf$kids_d6 + dataFLf$yr_96 + dataFLf$yr_97 ,bws=bwPf,regtype="lc", ckertype="epanechnikov"))
  pZm =  predict(npreg(formula = as.numeric(dataFLm$work) ~ dataFLm$benefit + dataFLm$married  + dataFLm$ed17 + dataFLm$ed18 + dataFLm$age  + dataFLm$reg_d1 + dataFLm$reg_d2 + dataFLm$reg_d3 + dataFLm$reg_d4 + dataFLm$reg_d5 + dataFLm$reg_d6 + dataFLm$reg_d7 + dataFLm$reg_d8 + dataFLm$reg_d9 + dataFLm$reg_d10 + dataFLm$reg_d11  + dataFLm$kids_d1 + dataFLm$kids_d2 + dataFLm$kids_d3 + dataFLm$kids_d4 + dataFLm$kids_d5 + dataFLm$kids_d6 + dataFLm$yr_96 + dataFLm$yr_97 ,bws=bwPm,regtype="lc", ckertype="epanechnikov"))
  


# Estimate parametric (probit) propensity score:

pZf_probit <- predict(glm(formula = dataFLf$work ~ dataFLf$benefit + dataFLf$married  + dataFLf$ed17 + dataFLf$ed18 + dataFLf$age  + dataFLf$reg_d1 + dataFLf$reg_d2 + dataFLf$reg_d3 + dataFLf$reg_d4 + dataFLf$reg_d5 + dataFLf$reg_d6 + dataFLf$reg_d7 + dataFLf$reg_d8 + dataFLf$reg_d9 + dataFLf$reg_d10 + dataFLf$reg_d11  + dataFLf$kids_d1 + dataFLf$kids_d2 + dataFLf$kids_d3 + dataFLf$kids_d4 + dataFLf$kids_d5 + dataFLf$kids_d6 + dataFLf$yr_96 + dataFLf$yr_97, family = binomial(link = "probit")), type="response")
pZm_probit <- predict(glm(formula=dataFLm$work ~ dataFLm$benefit + dataFLm$married  + dataFLm$ed17 + dataFLm$ed18 + dataFLm$age  + dataFLm$reg_d1 + dataFLm$reg_d2 + dataFLm$reg_d3 + dataFLm$reg_d4 + dataFLm$reg_d5 + dataFLm$reg_d6 + dataFLm$reg_d7 + dataFLm$reg_d8 + dataFLm$reg_d9 + dataFLm$reg_d10 + dataFLm$reg_d11  + dataFLm$kids_d1 + dataFLm$kids_d2 + dataFLm$kids_d3 + dataFLm$kids_d4 + dataFLm$kids_d5 + dataFLm$kids_d6 + dataFLm$yr_96 + dataFLm$yr_97, family = binomial(link = "probit")), type="response")

rm(data,data1,data2,datatmp)
str(dataFLf)
str(dataFLm)

#pZ = predict(glm(formula = dataFL$work ~ dataFL$benefit + factor(dataFL$married) + factor(dataFL$gender) + factor(dataFL$ed17) + factor(dataFL$ed18) + dataFL$age  + factor(dataFL$north) + factor(dataFL$midl) + factor(dataFL$east) + factor(dataFL$ldon) + factor(dataFL$south) + factor(dataFL$wales) + factor(dataFL$scot)  + factor(dataFL$kids_d1) + factor(dataFL$kids_d2) + factor(dataFL$kids_d3) + factor(dataFL$kids_d4) + factor(dataFL$kids_d5) + factor(dataFL$kids_d6) + factor(dataFL$yr_96) + factor(dataFL$yr_97) + factor(dataFL$yr_98) + factor(dataFL$yr_99) + factor(dataFL$yr_00),family=binomial(probit)), type="response")

dataFLf = data.frame(dataFLf,pZf,pZf_probit)
dataFLm = data.frame(dataFLm,pZm,pZm_probit)

# Selected Sample
dataFLfs =dataFLf[dataFLf$work==1,]
dataFLms =dataFLm[dataFLm$work==1,]

Yfs <- dataFLfs[,1]
Xfs <- dataFLfs[,4:(dim(dataFLfs)[2]-3)]
Pfs <- dataFLfs[,(dim(dataFLfs)[2]-1)]

Yms <- dataFLms[,1]
Xms <- dataFLms[,4:(dim(dataFLms)[2]-3)]
Pms <- dataFLms[,(dim(dataFLms)[2]-1)]



qfhat <- vector(mode="list",dim(as.matrix(tau_grid))[1])
qmhat <- vector(mode="list",dim(as.matrix(tau_grid))[1])


ehatf <- vector(mode="list",dim(as.matrix(tau_grid))[1])
ehatm <- vector(mode="list",dim(as.matrix(tau_grid))[1])

for (tau in 1:length(tau_grid)) {
  
  
  qqf <- matrix(0,length(Xfs),1,byrow=TRUE)
  qqm <- matrix(0,length(Xms),1,byrow=TRUE)
  
  
  
  for(i in 1:dim(Xfs)[1]) {
    Xfweight <- np::npksum(txdat=Xfs[,2:dim(Xfs)[2]],exdat=Xfs[i,2:dim(Xfs)[2]], leave.one.out=FALSE, bws=bwf , ckertype="epanechnikov", return.kernel.weights=TRUE)$kw
    qqf[i] <- rq(Yfs~1, weights=Xfweight, tau=tau_grid[[tau]], ci=FALSE)$coef[1.]
    
  }
  
  ehatf[[tau]] = Yfs-qqf
  qfhat[[tau]] <- qqf
  
  for(i in 1:dim(Xms)[1]) {
    
    Xmweight <- np::npksum(txdat=Xms[,2:dim(Xms)[2]],exdat=Xms[i,2:dim(Xms)[2]], leave.one.out=FALSE, bws=bwm , ckertype="epanechnikov", return.kernel.weights=TRUE)$kw
    qqm[i] <- rq(Yms~1, weights=Xmweight, tau=tau_grid[[tau]], ci=FALSE)$coef[1.]
    
    
  }
  
  ehatm[[tau]] = Yms-qqm
  qmhat[[tau]] <- qqm
  
}

## Estimation of Conditional (Unconditional) Distribution Function of quantile residuals:


Chatpf <- vector(mode="list",dim(as.matrix(tau_grid))[1])
Chatpm <- vector(mode="list",dim(as.matrix(tau_grid))[1])

for (tau in 1:length(tau_grid)) {
  
  BWef <- npregbw(formula=  Pfs ~ ehatf[[tau]],regtype="lc", ckertype="epanechnikov")$bw
  
  for (i in 1:length(Pfs)) {
    
    fphatf <- c(np::npksum(txdat=data.frame(Xfs[,2:dim(Xfs)[2]],ehatf[[tau]]),exdat=data.frame(Xfs[i,2:dim(Xfs)[2]],0), leave.one.out=FALSE, bandwidth.divide=TRUE, bws=c(bwCFf,BWef), ckertype="epanechnikov")$ksum)
    Chatpf[[tau]][i] <- (c(np::npksum(tydat=as.numeric(Pfs<= Pfs[i]),txdat=data.frame(Xfs[,2:dim(Xfs)[2]],ehatf[[tau]]),exdat=data.frame(Xfs[i,2:dim(Xfs)[2]],0), leave.one.out=FALSE, bandwidth.divide=TRUE, bws=c(bwCFf,BWef), ckertype="epanechnikov")$ksum))/fphatf
  }
  BWem <- npregbw(formula=  Pms ~ as.numeric(ehatm[[tau]]),regtype="lc", ckertype="epanechnikov")$bw
  
  for (i in 1:length(Pms)) {
    
    fphatm <- c(np::npksum(txdat=data.frame(Xms[,2:dim(Xms)[2]],ehatm[[tau]]),exdat=data.frame(Xms[i,2:dim(Xms)[2]],0), leave.one.out=FALSE, bandwidth.divide=TRUE, bws=c(bwCFm,BWem), ckertype="epanechnikov")$ksum)
    Chatpm[[tau]][i] <-  (c(np::npksum(tydat=as.numeric(Pms<=Pms[i]),txdat=data.frame(Xms[,2:dim(Xms)[2]],ehatm[[tau]]),exdat=data.frame(Xms[i,2:dim(Xms)[2]],0), leave.one.out=FALSE, bandwidth.divide=TRUE, bws=c(bwCFm,BWem), ckertype="epanechnikov")$ksum))/fphatm
  }
  
}

rm(ehatf,fphatf,ehatm,fphatm)


save.image("SelectedDataTmp95-97_NP.RData",compress="xz")


###################
## Execute Tests ##
###################





## FIRST TEST


IndXf <- IndMatX(Xfs[,2:dim(Xfs)[2]])
IndZf <- IndMatZ(Pfs)

T1IndZf <- vector(mode="list",dim(as.matrix(tau_grid))[1])

for (tau in 1:length(tau_grid)) {
  T1IndZf[[tau]] <- T1IndMatZ(Pfs,Chatpf[[tau]])
}

# Reset seed:
set.seed(12345)

# compute test (females):
CQ_stat_1f <- VBDNTest(Yfs,IndXf,IndZf,T1IndZf,qfhat, tau_grid, size=c(.1,.05), B=1000)

IndXm <- IndMatX(Xms[,2:dim(Xms)[2]])
IndZm <- IndMatZ(Pms)

T1IndZm <- vector(mode="list",dim(as.matrix(tau_grid))[1])

for (tau in 1:length(tau_grid)) {
  T1IndZm[[tau]] <- T1IndMatZ(Pms,Chatpm[[tau]])
}

# Reset seed:
set.seed(12345)

# compute test (males):
CQ_stat_1m <- VBDNTest(Yms,IndXm,IndZm,T1IndZm,qmhat, tau_grid, size=c(.1,.05), B=1000)

## SECOND TEST

## Kernel Function (2nd order Epanechnikov)
Kweights <- function(arg) {
  arg <- as.matrix(arg)
  return((3/4)*(1-(arg^2))*(as.numeric(abs(arg)<=1)))
}

# Reset seed:
set.seed(12345)

## Set the trimming and the bw parameter:
delta <- .98
hp    <- .02

# females: 
# compute kernel weights
#kweightsf <-  Kweights((Pfs-delta)/hp)

# compute test (females):
#CQ_stat_2f <- LocVBDNTest(Yfs,qfhat, kweightsf, size=c(.1,.05), tau_grid, B=1000)

#males
# compute kernel weights
kweightsm <- Kweights((Pms-delta)/hp)

# compute test:
CQ_stat_2m_1 <- LocVBDNTest(Yms,qmhat, kweightsm, size=c(.1,.05), tau_grid, B=1000)

## Set the trimming and the bw parameter:
delta <- .99
hp    <- .02

# compute kernel weights
kweightsm <- Kweights((Pms-delta)/hp)

# compute test (males):
CQ_stat_2m_2 <- LocVBDNTest(Yms,qmhat, kweightsm, size=c(.1,.05), tau_grid, B=1000)

## Set the trimming and the bw parameter:
delta <- 1
hp    <- .02

# compute kernel weights
kweightsm <- Kweights((Pms-delta)/hp)

# compute test (males):
CQ_stat_2m_3 <- LocVBDNTest(Yms,qmhat, kweightsm, size=c(.1,.05), tau_grid, B=1000)

## Set the trimming and the bw parameter:
delta <- 1
hp    <- .01

# compute kernel weights
kweightsm <- Kweights((Pms-delta)/hp)

# compute test (males):
CQ_stat_2m_4 <- LocVBDNTest(Yms,qmhat, kweightsm, size=c(.1,.05), tau_grid, B=1000)



## SAVE WORK SPACE
save.image("SelectedDataFINAL95-97_NP.RData")

# Total Running Time:
cat(paste("Running Time: ",floor(as.matrix(proc.time()-time)[3,]/60),"mins ",round((as.matrix(proc.time()-time)[3,]-60*floor(as.matrix(proc.time()-time)[3,]/60))),"sec",sep=""))

