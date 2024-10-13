setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./HAL/basis.fun.R")
source("./HAL/fit.hal.R")    
source("./Simulation/semi_markov_sim.R")
library("data.table")
library("Matrix")
library("glmnet")
library("dplyr")
library("stats")
library("pracma")
library("incase")
library("AalenJohansen")
library("nprobust")

##Initialize parameters
tstop <- 5 #terminal time of trajectory
nCut <- 13 #number of covariate partitions for HAL estimators
maxitTrapz <- 6
epsilon <- 10^(-7)
h <- 1/10
time.grid <- seq(from=0,to=tstop,by=h)
x.grid <- seq(from=-4,to=4,by=h)
len.x.grid <- length(x.grid)
nstepV <- tstop/h
nstepVSmall <- nstepV/5
newdata <- expand.grid(X=x.grid, tnNext=time.grid, tn=time.grid) %>% data.table()
hUnits <-  function(x){return(floor(x/h))} #note: will lead to round-off errors, but anyway have round-off error due to discretization
whichIndex <- function(s,t,x){1+ #first element
    hUnits(4+x)+ #how many steps to reach the desired value of X
    length(x.grid)*hUnits(t)+ #how many steps to reach the desired value of tnNext; for each increase of tnNext, one has to pass all the X-values
    length(x.grid)*length(time.grid)*hUnits(s) #how many steps to reach the desired value of tn; for each increase of tn, one has to pass all X and tnNext values
}
body(trapzfun)[[5]] <- "if (a == b) return(list(value = 0, iter = 0, error = 0))"

###Differential equation solver for estimand
f <- NULL

rk4 <- function(df, a, b, f0, n) {
  
  h = (b-a)/n
  
  f[1] = f0
  for (i in 1:n) {
    s = a + h * (i-1)
    k1 = h * df(s,f[i])
    k2 = h * df(s+0.5*h, f[i]+0.5*k1)
    k3 = h * df(s+0.5*h, f[i]+0.5*k2)
    k4 = h * df(s+h, f[i]+k3)
    
    f[i+1] = f[i]+1/6*(k1+2*k2+2*k3+k4)
  }
  return(f)
}

diffVi <- function(s,t,v,haz23,x){
  d1 = haz23(s,t,x)*v[1]-1
  return(c(d1))
}

Vi <- function(t0,s,haz23,x){
  V <- rk4(function(t, v) diffVi(t,v,haz23=haz23,s=s,x=x), a=tstop, b=t0+epsilon, f0=0, nstepVSmall)
  return(tail(V,n=1))
}

diffVa <- function(t,v,haz12,haz13,haz23, x){
  d1 = (haz12(0,t,x)+haz13(0,t,x))*v[1]-(t<=tstop)*haz12(0,t,x)*Vi(t,t,haz23,x)
  return(c(d1))
}

Va <- function(t0,haz12,haz13,haz23,x){
  V <- rk4(function(t,v) diffVa(t,v, haz12=haz12, haz13=haz13, haz23=haz23, x=x), a=tstop, b=t0+epsilon, f0=0, nstepV)
  return(tail(V,n=1))
}

VaVec <- function(t0,haz12,haz13,haz23,x){
  V <- rk4(function(t,v) diffVa(t,v, haz12=haz12, haz13=haz13, haz23=haz23, x=x), a=tstop, b=t0+epsilon, f0=0, nstepV)
  return(V)
}

haz12 <- function(s,t,x){transition_rate(1,2,s,t,x)}
haz13 <- function(s,t,x){transition_rate(1,3,s,t,x)}
haz23 <- function(s,t,x){transition_rate(2,3,s,t,x)}
haz14 <- function(s,t,x){transition_rate(1,4,s,t,x)}

##Compute Va oracle values over a grid
hVa <- (tstop-0)/nstepV
Va.time.grid <- tstop-hVa*(0:nstepV)
Va.value.grid <- expand.grid(t=Va.time.grid, X=x.grid)
Va.value.grid$VaTrue <- NA

for(m in 1:length(x.grid)){
  Va.value.grid$VaTrue[(1+(m-1)*length(Va.time.grid)):(m*length(Va.time.grid))] <- VaVec(0,haz12,haz13,haz23,x.grid[m])
}

whichIndexVa <- function(t,x){
  1+ #first element
  hUnits(tstop-t)+ #how many steps to reach the desired value of t
  length(time.grid)*hUnits(4+x) #how many steps to reach the desired value of x
}

###Helper functions

##Compute log-likelihood for censoring hazard in misspecified parametric model
misspecLoglike14 <- function(df,beta){
  
  occTerm <- df %>% 
    filter(yn==1 & ynNext==4) %>% 
    mutate(loghazard = beta[1]+beta[2]*tnNext+beta[3]*X) %>%
    summarise(loghazard=sum(loghazard)) %>%
    as.numeric()
  
  expoTerm <- df %>%
    filter(yn==1) %>%
    mutate(expo=1/beta[2]*exp(beta[1]+beta[2]*tnNext+beta[3]*X)-1/beta[2]*exp(beta[1]+beta[2]*0+beta[3]*X)  ) %>% 
    summarise(expo = sum(expo)) %>%
    as.numeric()
  
  return(occTerm-expoTerm)
}


##Compute pseudo-outcomes for Df2 based on estimators based on Df1  
pseudoOutcomes <- function(Df1,Df2,seed){
  
  ##Fit HAL
  fit12 <- fit.hal(fixed.covars=c("X","tnNext"), time.covars=c(), Df1[yn == 1,], V=5, cut=nCut, seed=seed,
                   delta.var="ynNext", delta.value=2, time.var="tnNext", browse1=FALSE)
  
  gc()
  
  fit13 <- fit.hal(fixed.covars=c("X","tnNext"), time.covars=c(), Df1[yn == 1,], V=5, cut=nCut, seed=seed,
                   delta.var="ynNext", delta.value=3, time.var="tnNext", browse1=FALSE)
  
  gc()
  
  fit14 <- fit.hal(fixed.covars=c("X","tnNext"), time.covars=c(), Df1[yn == 1,], V=5, cut=nCut, seed=seed,
                   delta.var="ynNext", delta.value=4, time.var="tnNext", browse1=FALSE)
  
  gc()
  
  fit23 <- fit.hal(fixed.covars=c("X","tnNext","tn"), time.covars=c(), Df1[yn == 2,], V=5, cut=nCut, seed=seed,
                   delta.var="ynNext", delta.value=3, time.var="tnNext", browse1=FALSE)
  
  gc()
  
  fitHAL <-  function(i, j){
    if(i==1 && j==2){fit12}
    else if(i==1 && j==3){fit13}
    else if(i==1 && j==4){fit14}
    else if(i==2 && j==3){fit23}
    else{error("error: not a valid transition.")}
  }
  
  predHALDf <- function(i,j,newdata){
    fit <- fitHAL(i,j)
    coef <- coef(fit$fit, s = fit$cve$lambda.min)
    newdataMatrix <- Matrix(model.matrix(fit$model.formula.X, data=newdata), sparse=FALSE)
    as.numeric(exp(coef[1]+newdataMatrix %*% coef[-c(1)]))
  }
  
  hazHALMatrix <- newdata
  hazHALMatrix$haz12 <- predHALDf(1,2,newdata)
  hazHALMatrix$haz13 <- predHALDf(1,3,newdata)
  hazHALMatrix$haz14 <- predHALDf(1,4,newdata)
  hazHALMatrix$haz23 <- predHALDf(2,3,newdata)
  
  haz12HAL <- function(s,t,x){hazHALMatrix$haz12[whichIndex(s,t,x)]}
  haz13HAL <- function(s,t,x){hazHALMatrix$haz13[whichIndex(s,t,x)]}
  haz14HAL <- function(s,t,x){hazHALMatrix$haz14[whichIndex(s,t,x)]}
  haz23HAL <- function(s,t,x){hazHALMatrix$haz23[whichIndex(s,t,x)]}
  haz24HAL <- function(s,t,x){0}
  surv14HAL <- function(t1,t2,s,x){exp(-trapzfun(f=function(t){haz14HAL(s,t,x)}, t1, t2, maxit=maxitTrapz)$value)}
  surv24HAL <- function(t1,t2,s,x){1}
  
  ##Fit misspecified censoring hazard
  misspecParam14 <- optim(c(-0.1,-0.1,-0.1),fn=misspecLoglike14, df=dt, control = list(fnscale = -1))$par
  haz14Misspec <- function(s,t,x){exp(misspecParam14[1]+misspecParam14[2]*t+misspecParam14[3]*x)}
  haz24Misspec <- function(s,t,x){0}
  surv14Misspec <- function(t1,t2,s,x){exp(-(1/misspecParam14[2]*exp(misspecParam14[1]+misspecParam14[2]*t2+misspecParam14[3]*x)-1/misspecParam14[2]*exp(misspecParam14[1]+misspecParam14[2]*t1+misspecParam14[3]*x)))}
  surv24Misspec <- function(t1,t2,s,x){1}
  
  ##Oracle censoring 
  surv14Oracle <- function(t1,t2,s,x){exp(-(t2-t1)*transition_rate(1,4,0,0,x))}
  surv24Oracle <- function(t1,t2,s,x){1}
  
  ##Compute Va HAL plug-in over a grid
  for(m in 1:length(x.grid)){
    Va.value.grid$VaHAL[(1+(m-1)*length(Va.time.grid)):(m*length(Va.time.grid))] <- VaVec(0,haz12HAL,haz13HAL,haz23HAL,x.grid[m])
  }
  
  ##Output
  Df2Pseudo <- Df2 %>% group_by(id) %>% mutate(
    IPCWTerm = in_case(
      max(yn)==2 & max(ynNext)!=4 ~ (max(tnNext)-max(tn))/(surv14HAL(0,max(tn),0,max(X))*surv24HAL(max(tn),max(tnNext),max(tn),max(X))),
      TRUE ~ 0), 
    IPCWTermMisspec = in_case(
      max(yn)==2 & max(ynNext)!=4 ~ (max(tnNext)-max(tn))/(surv14Misspec(0,max(tn),0,max(X))*surv24Misspec(max(tn),max(tnNext),max(tn),max(X))),
      TRUE ~ 0),
    IPCWTermOracle = in_case(
      max(yn)==2 & max(ynNext)!=4 ~ (max(tnNext)-max(tn))/(surv14Oracle(0,max(tn),0,max(X))*surv24Oracle(max(tn),max(tnNext),max(tn),max(X))),
      TRUE ~ 0)) %>% mutate(
    JumpTerm = in_case(
      max(yn)==1 & max(ynNext)==4 ~ Va.value.grid$VaHAL[whichIndexVa(max(tnNext),max(X))]/surv14HAL(0,max(tnNext),0,max(X)),
      max(yn)==2 & max(ynNext)==4 ~ (max(tnNext)-max(tn)+Vi(max(tnNext),max(tn),haz23HAL,max(X)))/(surv14HAL(0,max(tn),0,max(X))*surv24HAL(max(tn),max(tnNext),max(tn),max(X))),
      TRUE ~ 0),
    JumpTermMisspec = in_case(
      max(yn)==1 & max(ynNext)==4 ~ Va.value.grid$VaHAL[whichIndexVa(max(tnNext),max(X))]/surv14Misspec(0,max(tnNext),0,max(X)),
      max(yn)==2 & max(ynNext)==4 ~ (max(tnNext)-max(tn)+Vi(max(tnNext),max(tn),haz23HAL,max(X)))/(surv14Misspec(0,max(tn),0,max(X))*surv24Misspec(max(tn),max(tnNext),max(tn),max(X))),
      TRUE ~ 0),
    JumpTermOracle = in_case(
      max(yn)==1 & max(ynNext)==4 ~ Va.value.grid$VaTrue[whichIndexVa(max(tnNext),max(X))]/surv14Oracle(0,max(tnNext),0,max(X)),
      max(yn)==2 & max(ynNext)==4 ~ (max(tnNext)-max(tn)+Vi(max(tnNext),max(tn),haz23,max(X)))/(surv14Oracle(0,max(tn),0,max(X))*surv24Oracle(max(tn),max(tnNext),max(tn),max(X))),
      TRUE ~ 0)) %>% mutate(
    CompensatorTerm = in_case(
      max(yn)==1 ~ trapzfun(f=function(t){Va.value.grid$VaHAL[whichIndexVa(t,max(X))]*haz14HAL(0,t,max(X))/surv14HAL(0,t,0,max(X))}, 0, max(tnNext), maxit=maxitTrapz)$value,
      max(yn)==2 ~ trapzfun(f=function(t){Va.value.grid$VaHAL[whichIndexVa(t,max(X))]*haz14HAL(0,t,max(X))/surv14HAL(0,t,0,max(X))}, 0, max(tn), maxit=maxitTrapz)$value, 
      TRUE ~ "ERROR"),
    CompensatorTermMisspec = in_case(
      max(yn)==1 ~ trapzfun(f=function(t){Va.value.grid$VaHAL[whichIndexVa(t,max(X))]*haz14Misspec(0,t,max(X))/surv14Misspec(0,t,0,max(X))}, 0, max(tnNext), maxit=maxitTrapz)$value,
      max(yn)==2 ~ trapzfun(f=function(t){Va.value.grid$VaHAL[whichIndexVa(t,max(X))]*haz14Misspec(0,t,max(X))/surv14Misspec(0,t,0,max(X))}, 0, max(tn), maxit=maxitTrapz)$value, 
      TRUE ~ "ERROR"),
    CompensatorTermOracle = in_case(
      max(yn)==1 ~ trapzfun(f=function(t){Va.value.grid$VaTrue[whichIndexVa(t,max(X))]*haz14(0,t,max(X))/surv14Oracle(0,t,0,max(X))}, 0, max(tnNext), maxit=maxitTrapz)$value,
      max(yn)==2 ~ trapzfun(f=function(t){Va.value.grid$VaTrue[whichIndexVa(t,max(X))]*haz14(0,t,max(X))/surv14Oracle(0,t,0,max(X))}, 0, max(tn), maxit=maxitTrapz)$value, 
      TRUE ~ "ERROR")) %>% mutate(
    CompensatorTerm=as.numeric(CompensatorTerm), 
    CompensatorTermMisspec=as.numeric(CompensatorTermMisspec),
    CompensatorTermOracle=as.numeric(CompensatorTermOracle)) %>% mutate(
    Ystar = IPCWTerm+JumpTerm-CompensatorTerm, 
    YstarMisspec = IPCWTermMisspec+JumpTermMisspec-CompensatorTermMisspec,
    YstarOracle = IPCWTermOracle+JumpTermOracle-CompensatorTermOracle)
  
  gc()
  
  return(Df2Pseudo)
}


###Simulate
seeds <- 1:500
nsims <- c(1000,5000,10000,30000)

namesRes <- c("seed","n","metric","estimator","value")
namesmetrics <- c("RMISE","MIAE","MD","SqErr0","AbsErr0","SqErrNeg1","AbsErrorNeg1")
namesEstimators <- c("CAJ","HAL","DR","DR.Misspecified", "DR.Oracle", "IPCW", "IPCW.Misspecified", "IPCW.Oracle")
indexRes <- 1
dfRes <- data.frame(matrix(ncol = length(namesRes), nrow = length(seeds)*length(nsims)*length(namesmetrics)*length(namesEstimators)))
colnames(dfRes) <- namesRes

namesResCI <- c("seed","n","X","estimator","Pred","SE")
namesEstimatorsCI <- c("DR","DR.Misspecified", "DR.Oracle", "IPCW", "IPCW.Misspecified", "IPCW.Oracle")
indexResCI <- 1
dfResCI <- data.frame(matrix(ncol = length(namesResCI), nrow = length(seeds)*length(nsims)*length(x.grid)*length(namesEstimatorsCI)))
colnames(dfResCI) <- namesResCI

for(seed in seeds){
for(nsim in nsims){

set.seed(seed)
X <- runif(nsim,min=-4,max=4) #simulate X

set.seed(seed)
simList <- simulator_n(1, 0, 0, tstop, X, TRUE)

set.seed(seed)
dataRaw <- simulator_n(1, 0, 0, tstop, X, FALSE)

data <- dataRaw %>% filter(!is.na(tn))  #remove "all-NA" rows
dt <- data %>%
  filter(!(yn %in% c(3,4))) %>% #no transitions possible after transition to state 3 or 4
  mutate(tnNext=ifelse(is.na(tnNext),tstop,tnNext),
         ynNext=ifelse(is.na(ynNext),5,ynNext) ) %>% #add transition to state 5 as administrative censoring (state 4 is drop-out censoring)
  data.table()


###Compute estimators

##Compute true value
vaPred <- Va.value.grid %>% filter(t==0) %>% select(VaTrue) %>% unlist() %>% as.numeric()

##Compute Conditional Aalen-Johansen plug-in estimator
vaAJPred <- sapply(x.grid, FUN=function(x){
  AJfit <- aalen_johansen(simList, x = x, alpha=0.2)
  AJtime <- AJfit$t
  AJprob <- unlist(lapply(AJfit$p, FUN = function(L) L[2]))
  sum(diff(AJtime)*head(AJprob,-1))
})

gc()

##Fit HAL
fit12Full <- fit.hal(fixed.covars=c("X","tnNext"), time.covars=c(), dt[yn == 1,], V=5, cut=nCut, seed=seed,
                     delta.var="ynNext", delta.value=2, time.var="tnNext", browse1=FALSE)

gc()

fit13Full <- fit.hal(fixed.covars=c("X","tnNext"), time.covars=c(), dt[yn == 1,], V=5, cut=nCut, seed=seed,
                     delta.var="ynNext", delta.value=3, time.var="tnNext", browse1=FALSE)

gc()

fit23Full <- fit.hal(fixed.covars=c("X","tnNext","tn"), time.covars=c(), dt[yn == 2,], V=5, cut=nCut, seed=seed,
                     delta.var="ynNext", delta.value=3, time.var="tnNext", browse1=FALSE)

gc()

fitHALFull <-  function(i, j){
  if(i==1 && j==2){fit12Full}
  else if(i==1 && j==3){fit13Full}
  else if(i==2 && j==3){fit23Full}
  else{error("error: not a valid transition.")}
}

predHALFullDf <- function(i,j,newdata){
  fit <- fitHALFull(i,j)
  coef <- coef(fit$fit, s = fit$cve$lambda.min)
  newdataMatrix <- Matrix(model.matrix(fit$model.formula.X, data=newdata), sparse=FALSE)
  as.numeric(exp(coef[1]+newdataMatrix %*% coef[-c(1)]))
}

##Calculate HAL hazard values over a grid
hazHALFullMatrix <- newdata
hazHALFullMatrix$haz12 <- predHALFullDf(1,2,newdata)
hazHALFullMatrix$haz13 <- predHALFullDf(1,3,newdata)
hazHALFullMatrix$haz23 <- predHALFullDf(2,3,newdata)

haz12HALFull <- function(s,t,x){hazHALFullMatrix$haz12[whichIndex(s,t,x)]}
haz13HALFull <- function(s,t,x){hazHALFullMatrix$haz13[whichIndex(s,t,x)]}
haz23HALFull <- function(s,t,x){hazHALFullMatrix$haz23[whichIndex(s,t,x)]}

vaHALFullPred <- sapply(x.grid, FUN=function(x){Va(0,haz12HALFull,haz13HALFull,haz23HALFull,x)})

##Debiased plug-in HAL cross-fitted with local polynomial regression
percentD1 <- 0.5
Df1 <- dt[which(dt$id <= nsim*percentD1),]
Df2 <- dt[which(dt$id > nsim*percentD1),]

Df2Pseudo <- pseudoOutcomes(Df1,Df2,seed)
Df1Pseudo <- pseudoOutcomes(Df2,Df1,seed)


Df2PseudoInput <- Df2Pseudo %>% group_by(id) %>% summarise(X=max(X), IPCWTerm = max(IPCWTerm), Ystar=max(Ystar), 
                                                           IPCWTermMisspec=max(IPCWTermMisspec), YstarMisspec=max(YstarMisspec),
                                                           IPCWTermOracle=max(IPCWTermOracle), YstarOracle=max(YstarOracle))
Df1PseudoInput <- Df1Pseudo %>% group_by(id) %>% summarise(X=max(X), IPCWTerm = max(IPCWTerm), Ystar=max(Ystar), 
                                                           IPCWTermMisspec=max(IPCWTermMisspec), YstarMisspec=max(YstarMisspec),
                                                           IPCWTermOracle=max(IPCWTermOracle), YstarOracle=max(YstarOracle))

#DR
vaDRDf2Pred <- lprobust(Df2PseudoInput$Ystar,Df2PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
vaDRDf1Pred <- lprobust(Df1PseudoInput$Ystar,Df1PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
gc()

#DR misspecified
vaDRDf2MisspecPred <- lprobust(Df2PseudoInput$YstarMisspec,Df2PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
vaDRDf1MisspecPred <- lprobust(Df1PseudoInput$YstarMisspec,Df1PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
gc()

#DR oracle
vaDRDf2OraclePred <- lprobust(Df2PseudoInput$YstarOracle,Df2PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
vaDRDf1OraclePred <- lprobust(Df1PseudoInput$YstarOracle,Df1PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
gc()

#IPCW
vaIPCWDf2Pred <- lprobust(Df2PseudoInput$IPCWTerm,Df2PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
vaIPCWDf1Pred <- lprobust(Df1PseudoInput$IPCWTerm,Df1PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
gc()

#IPCW misspecified
vaIPCWDf2MisspecPred <- lprobust(Df2PseudoInput$IPCWTermMisspec,Df2PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
vaIPCWDf1MisspecPred <- lprobust(Df1PseudoInput$IPCWTermMisspec,Df1PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
gc()

#IPCW oracle
vaIPCWDf2OraclePred <- lprobust(Df2PseudoInput$IPCWTermOracle,Df2PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
vaIPCWDf1OraclePred <- lprobust(Df1PseudoInput$IPCWTermOracle,Df1PseudoInput$X,eval=x.grid,h=(nsim/5000)^(-1/4.5)*1.3)
gc()

##Save results estimators
save.est <- function(name,f){
  dfRes[indexRes,] <<- c(seed,nsim,name,"CAJ",f(vaAJPred)) 
  indexRes <<- indexRes+1
  
  dfRes[indexRes,] <<- c(seed,nsim,name,"HAL",f(vaHALFullPred)) 
  indexRes <<- indexRes+1
  
  dfRes[indexRes,] <<- c(seed,nsim,name,"DR.HAL",f(vaDRDf1Pred$Estimate[,"tau.us"]*0.5+vaDRDf2Pred$Estimate[,"tau.us"]*0.5)) 
  indexRes <<- indexRes+1
  
  dfRes[indexRes,] <<- c(seed,nsim,name,"DR.Missspecified",f(vaDRDf1MisspecPred$Estimate[,"tau.us"]*0.5+vaDRDf2MisspecPred$Estimate[,"tau.us"]*0.5)) 
  indexRes <<- indexRes+1
  
  dfRes[indexRes,] <<- c(seed,nsim,name,"DR.Oracle",f(vaDRDf1OraclePred$Estimate[,"tau.us"]*0.5+vaDRDf2OraclePred$Estimate[,"tau.us"]*0.5)) 
  indexRes <<- indexRes+1
  
  dfRes[indexRes,] <<- c(seed,nsim,name,"IPCW.HAL",f(vaIPCWDf1Pred$Estimate[,"tau.us"]*0.5+vaIPCWDf2Pred$Estimate[,"tau.us"]*0.5)) 
  indexRes <<- indexRes+1
  
  dfRes[indexRes,] <<- c(seed,nsim,name,"IPCW.Missspecified",f(vaIPCWDf1MisspecPred$Estimate[,"tau.us"]*0.5+vaIPCWDf2MisspecPred$Estimate[,"tau.us"]*0.5)) 
  indexRes <<- indexRes+1
  
  dfRes[indexRes,] <<- c(seed,nsim,name,"IPCW.Oracle",f(vaIPCWDf1OraclePred$Estimate[,"tau.us"]*0.5+vaIPCWDf2OraclePred$Estimate[,"tau.us"]*0.5)) 
  indexRes <<- indexRes+1
  
  gc()
}

#root mean integrated squared error (RMISE)
rmise <- function(x){sqrt(sum(diff(x.grid)*head((vaPred-x)^2,-1),na.rm=TRUE))}
save.est("RMISE",rmise)

#Mean integrated absolute error (MIAE)
miae <- function(x){sum(abs(diff(x.grid)*head((vaPred-x),-1)),na.rm=TRUE)}
save.est("MIAE",miae)

#Maximum distance (MD)
md <- function(x){max(abs(vaPred-x), na.rm=TRUE)}
save.est("MD",md)

#Squared error at X=0 and X=-1
x0Idx <- which(x.grid==0)
xNeg1Idx <- which(x.grid==-1)

sqErr0 <- function(x){(vaPred[x0Idx]-x[x0Idx])^2}
sqErrNeg1 <- function(x){(vaPred[xNeg1Idx]-x[xNeg1Idx])^2}
save.est("SqErr0",sqErr0)
save.est("SqErrNeg1",sqErrNeg1)

#Absolute error at X=0 and X=-1
absErr0 <- function(x){abs(vaPred[x0Idx]-x[x0Idx])}
absErrNeg1 <- function(x){abs(vaPred[xNeg1Idx]-x[xNeg1Idx])}
save.est("AbsErr0",absErr0)
save.est("AbsErrNeg1",absErrNeg1)

##Save results CI
save.CI <- function(name,fit1,fit2){
  dfResCI[indexResCI:(indexResCI+len.x.grid-1),] <<- data.frame(seed=rep(seed,len.x.grid),
                                                               n=rep(nsim,len.x.grid),
                                                               X=x.grid,
                                                               estimator=rep(name,len.x.grid),
                                                               Pred=0.5*fit1$Estimate[,"tau.us"]+0.5*fit2$Estimate[,"tau.us"],
                                                               SE=0.5*fit1$Estimate[,"se.us"]+0.5*fit2$Estimate[,"se.us"]) 
  indexResCI <<- indexResCI+len.x.grid
  
  gc()
}
save.CI("DR",vaDRDf1Pred,vaDRDf2Pred)
save.CI("DR.Missspecified",vaDRDf1MisspecPred,vaDRDf2MisspecPred)
save.CI("DR.Oracle",vaDRDf1OraclePred,vaDRDf2OraclePred)
save.CI("IPCW",vaIPCWDf1Pred,vaIPCWDf2Pred)
save.CI("IPCW.Missspecified",vaIPCWDf1MisspecPred,vaIPCWDf2MisspecPred)
save.CI("IPCW.Oracle",vaIPCWDf1OraclePred,vaIPCWDf2OraclePred)

##To save partial results, uncomment here:
#if((seed == 1 | seed %% 5 == 0) & nsim == max(nsims)){
#  save(dfRes, file = paste0("./Results/dfRes_",seed,".Rda"))
#  save(dfResCI, file = paste0("./Results/dfResCI_",seed,".Rda"))
#}

gc()

if(seed == 1 & nsim == max(nsims)){
  #Estimand plot
  png("Figures/EstimandPlot.png", width = 6, height = 4, units = 'in', res = 400)
  par(mar = c(2.5, 2.5, 0.5, 0))
  par(cex.lab=0.9, cex.axis=0.9)
  par(mgp=c(1.7, 0.5, 0)) 
  
  plot(x.grid,vaPred,type="l", ylim=c(0.77,1.28), lwd=2, xlab="W", ylab="Expected illness duration")
  lines(x.grid,vaAJPred, col="burlywood2")
  lines(x.grid,vaHALFullPred, col="burlywood4")
  lines(x.grid,vaDRDf1Pred$Estimate[,"tau.us"]*0.5+vaDRDf2Pred$Estimate[,"tau.us"]*0.5, col="ForestGreen")
  lines(x.grid,vaDRDf1MisspecPred$Estimate[,"tau.us"]*0.5+vaDRDf2MisspecPred$Estimate[,"tau.us"]*0.5, col="ForestGreen", lty=2)
  lines(x.grid,vaDRDf1OraclePred$Estimate[,"tau.us"]*0.5+vaDRDf2OraclePred$Estimate[,"tau.us"]*0.5, col="ForestGreen", lty=3)
  lines(x.grid,vaIPCWDf1Pred$Estimate[,"tau.us"]*0.5+vaIPCWDf2Pred$Estimate[,"tau.us"]*0.5, col="steelblue1")
  lines(x.grid,vaIPCWDf1MisspecPred$Estimate[,"tau.us"]*0.5+vaIPCWDf2MisspecPred$Estimate[,"tau.us"]*0.5, col="steelblue1", lty=2)
  lines(x.grid,vaIPCWDf1OraclePred$Estimate[,"tau.us"]*0.5+vaIPCWDf2OraclePred$Estimate[,"tau.us"]*0.5, col="steelblue1", lty=3)
  legend("topright", 
         legend=c("True value", "Plug-in CAJ", "Plug-in HAL", "", "", "", "DR HAL", "DR misspecified", "DR oracle", "IPCW", "IPCW misspecified", "IPCW oracle"),
         col=c("black", "burlywood2", "burlywood4", "white", "white", "white", "ForestGreen", "ForestGreen", "ForestGreen", "steelblue1", "steelblue1", "steelblue1"), 
         lty=c(1,1,1,0,0,0,1,2,3,1,2,3),
         cex=0.9, 
         lwd=c(2,1,1,1,1,1,1,1,1,1,1,1),
         box.lty=0,
         inset=.001,
         ncol=2,
         bg="transparent")
  
  dev.off()

  ##Hazard plot
  png("Figures/HazardPlot.png", width = 6, height = 4, units = 'in', res = 400)
  par(mfrow=c(3,2))
  par(cex.lab=0.9, cex.axis=0.9)
  par(mar = c(2, 2, 1, 1))
  par(mgp=c(0,0.7,0))
  
  time.grid.plot <- seq(from=0,to=tstop-2,by=h)
  
  plot(x.grid,haz12HALFull(0,2,x.grid),ylim=c(0.23,0.46),type="l", col="blue", xlab="", ylab="")
  legend("topright", legend=expression(paste('μ'[12],'(2,w)')), bty = "n")
  lines(x.grid,transition_rate(1,2,0,2,x.grid)) 
  
  plot(time.grid,haz12HALFull(0,time.grid,-1),type="l", col="blue", ylim=c(0.28,0.4), xlab="", ylab="") #plot(time.grid,haz12HALFull(0,time.grid,0),type="l", col="blue", ylim=c(0.33,0.45), xlab="", ylab="")
  legend("topright", legend=expression(paste('μ'[12],'(t,-1)')), bty = "n")
  lines(time.grid,transition_rate(1,2,0,time.grid,-1)) 
  
  plot(x.grid,haz13HALFull(0,2,x.grid),ylim=c(0.08,0.16), type="l", col="blue", xlab="", ylab="")
  legend("topright", legend=expression(paste('μ'[13],'(2,w)')), bty = "n")
  lines(x.grid,transition_rate(1,3,0,2,x.grid)) 
  
  plot(x.grid,haz23HALFull(0.5,2,x.grid), type="l", col="blue", ylim=c(0.28,0.4), xlab="", ylab="")
  legend("topright", legend=expression(paste('μ'[23],'(2,0.5,w)')), bty = "n")
  lines(x.grid,lapply(x.grid, function(x)transition_rate(2,3,0.5,2,x))) 
  
  plot(time.grid.plot,haz23HALFull(2,2+time.grid.plot,-1), type="l", col="blue", ylim=c(0,1.1), xlab="", ylab="")
  legend("topright", legend=expression(paste('μ'[23],'(2,2+t,-1)')), bty = "n")
  lines(time.grid.plot,sapply(time.grid.plot, FUN=function(t) transition_rate(2,3,2,2+t,-1)) ) 
  
  fit14Full <- fit.hal(fixed.covars=c("X","tnNext"), time.covars=c(), dt[yn == 1,], V=5, cut=nCut, seed=seed,
                       delta.var="ynNext", delta.value=4, time.var="tnNext", browse1=FALSE)
  coef14 <- coef(fit14Full$fit, s = fit14Full$cve$lambda.min)
  newdataMatrix14 <- Matrix(model.matrix(fit14Full$model.formula.X, data=newdata %>% filter(tnNext==2 & tn==0)), sparse=FALSE)
  pred14 <- as.numeric(exp(coef14[1]+newdataMatrix14 %*% coef14[-c(1)]))
  plot(x.grid,pred14, type="l", col="blue", ylim=c(0.19,0.415), xlab="", ylab="")
  legend("topright", legend=expression(paste('1{Z(2)=1}',gamma,'(2,w)')), bty = "n")
  lines(x.grid,transition_rate(1,4,0,2,x.grid))
  
  dev.off()
  
  ##CI plot
  png("Figures/CIPlot.png", width = 8, height = 4, units = 'in', res = 400)
  par(mar = c(2.5, 2.5, 0, 0))
  par(cex.lab=1, cex.axis=1)
  par(mgp=c(1.7, 0.5, 0)) 
  
  plot(x.grid,vaPred,type="l", ylim=c(0.8,1.21), lwd=2, xlab="W", ylab="Expected illness duration")
  lines(x.grid,vaDRDf1Pred$Estimate[,"tau.us"]*0.5+vaDRDf2Pred$Estimate[,"tau.us"]*0.5, col="ForestGreen")
  SeConstant <- 1/sqrt(2)
  lines(x.grid,(vaDRDf1Pred$Estimate[,"tau.us"]+1.96*SeConstant*vaDRDf1Pred$Estimate[,"se.us"])*0.5+(vaDRDf2Pred$Estimate[,"tau.us"]+1.96*SeConstant*vaDRDf2Pred$Estimate[,"se.us"])*0.5, col="blue", lty=2)
  lines(x.grid,(vaDRDf1Pred$Estimate[,"tau.us"]-1.96*SeConstant*vaDRDf1Pred$Estimate[,"se.us"])*0.5+(vaDRDf2Pred$Estimate[,"tau.us"]-1.96*SeConstant*vaDRDf2Pred$Estimate[,"se.us"])*0.5, col="blue", lty=2)
  legend("topright", legend=c("True value","DR HAL","95% confidence band"),
         col=c("black", "ForestGreen", "blue"), 
         lty=c(1,1,2),
         cex=1, 
         lwd=c(2,1,1),
         box.lty=0,
         inset=.001,
         bg="transparent")
  
  dev.off()
}

}
}

save(dfRes, file = paste0("./Results/dfRes.Rda"))
save(dfResCI, file = paste0("./Results/dfResCI.Rda"))
