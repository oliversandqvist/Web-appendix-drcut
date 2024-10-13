setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./HAL/basis.fun.R")
source("./HAL/fit.hal.R")    
source("./Simulation/semi_markov_sim.R")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(Hmisc)

##Initialize parameters
tstop <- 5 #terminal time of trajectory
nCut <- 13 #number of covariate partitions for HAL estimators
maxitTrapz <- 6
epsilon <- 10^(-7)
h <- 1/10
time.grid <- seq(from=0,to=tstop,by=h)
x.grid <- seq(from=-4,to=4,by=h)
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

Va.value.grid <- Va.value.grid %>% filter(t==0) %>% select(X,VaTrue)

#predictive performance
load("./Results/dfRes.Rda")

dfRes <- dfRes %>% 
  filter(!is.na(seed)) %>%
  mutate(value=as.numeric(value)) %>%
  mutate(n=factor(n, levels=c('1000','5000','10000','30000'))) %>%
  mutate(metric=factor(metric,levels=c('RMISE','MIAE','MD','SqErr0','SqErrNeg1','AbsErr0','AbsErrNeg1')))

##Violinplot for RMISE
png("Figures/ViolinPlot.png", width = 12, height = 5.5, units = 'in', res = 600)
par(mar = c(0, 0, 0, 0))
par(mgp=c(0, 0, 0)) 

ggplot(aes(y = value, x = estimator),
       data = dfRes %>% filter(metric == "RMISE") %>% filter(value < 1)) +
  theme_light() +
  geom_violin(color="black") +
  facet_grid(. ~ n) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + 
  ylab(expression(paste(L^2,"([-4,4],",lambda,")"," error"))) +
  xlab("") +
  theme(text = element_text(size = 13),
        strip.text = element_text(size=13, face="bold")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black", fun.args = list(mult = 1))

dev.off()

#inference
load("./Results/dfResCI.Rda")

dfResCI <- dfResCI %>% left_join(Va.value.grid, by=c("X"))

##CI empirical coverages with DR HAL and DR Oracle using standard errors outputted by lprobust

SeConstant <- 1/sqrt(2)

dfResCIPlot <- dfResCI %>% 
  filter(estimator %in% c("DR","DR.Oracle") ) %>%
  group_by(X,n,estimator) %>%
  summarise(pct99 = 100-sum(abs(Pred-VaTrue) > SeConstant*2.58*SE)/n()*100, 
            pct95 = 100-sum(abs(Pred-VaTrue) > SeConstant*1.96*SE)/n()*100,
            pct90 = 100-sum(abs(Pred-VaTrue) > SeConstant*1.65*SE)/n()*100) %>%
  pivot_longer(cols=starts_with("pct"),names_to="significance")

png("Figures/CoveragePlot.png", width = 12, height = 7, units = 'in', res = 600)
par(mar = c(0, 0, 0, 0))
par(mgp=c(0, 0, 0)) 

ggplot(aes(y = value, x = X, color=significance, linetype=estimator),
       data = dfResCIPlot) +
  theme_light() +
  geom_line(lwd=1.25) +
  geom_hline(yintercept=99, color="gray80", lty=2, lwd=1.2) +
  geom_hline(yintercept=95, color="gray80", lty=2, lwd=1.2) +
  geom_hline(yintercept=90, color="gray80", lty=2, lwd=1.2) +
  facet_grid(. ~ n) +
  scale_linetype_manual(values=c(1,3)) +
  ylab(expression(paste("Empirical coverage"))) +
  xlab("W") +
  theme(text = element_text(size = 24),
        strip.text = element_text(size=22, face="bold")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom",
        legend.key.size = unit(1.25, 'cm'),
        legend.title = element_text(size=22),
        legend.text = element_text(size = 20)) +
  scale_color_discrete(name = "Significance", labels = c("90%", "95%", "99%")) +
  scale_linetype_discrete(name = "Estimator", labels = c("DR", "DR oracle")) 

dev.off()

##Histogram over predicted values at X=-1 for HAL and DR, and overlay gaussian distribution based on approximation for oracle outcomes

png("Figures/HistogramPlot.png", width = 10, height = 5, units = 'in', res = 600)
par(mar = c(2, 2, 1.5, 1))
par(mgp=c(3, 1, 0)) 
par(mfrow = c(1, 2))
par(cex.lab=1.3, cex.axis=1.3, cex.main=1.4)

predOracle <- dfResCI %>% filter(X==-1 & estimator == "DR.Oracle" & abs(Pred) < 3 & n==30000) %>% select(Pred) %>% unlist()
histDRDf <- dfResCI %>% filter(X==-1 & estimator == "DR" & abs(Pred) < 3 & n==30000) 
predDR <- histDRDf %>% select(Pred) %>% unlist()
predDRMiss <- dfResCI %>% filter(X==-1 & estimator == "DR.Missspecified" & abs(Pred) < 3 & n==30000) %>% select(Pred) %>% unlist()

hist(predDRMiss, breaks=30, probability=TRUE, xlim=c(0.9,1.15), ylim=c(0,21), xlab="", ylab="", main="DR misspecified IPCW, 500 simulations", col="steelblue2")
curve(dnorm(x, mean=mean(predOracle), sd=sd(predOracle)),
      col="red", lwd=2.5, add=TRUE, yaxt="n")
abline(v=Va.value.grid%>%filter(X==-1)%>%select(VaTrue)%>%unlist(), col="red", lwd=2, lty=2)

hist(predDR, breaks=30, probability=TRUE, xlim=c(0.9,1.15), ylim=c(0,21), xlab="", ylab="", main="DR, 500 simulations", col="steelblue2")
curve(dnorm(x, mean=mean(predOracle), sd=sd(predOracle)),
      col="red", lwd=2.5, add=TRUE, yaxt="n")
abline(v=Va.value.grid%>%filter(X==-1)%>%select(VaTrue)%>%unlist(), col="red", lwd=2, lty=2)

dev.off()

