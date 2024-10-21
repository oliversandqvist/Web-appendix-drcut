library(foreign)
library(tidyverse)
library(xgboost)
library(nprobust)
library(rdrobust)
library(rdmulti)
library(stringr)
library(gridExtra)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

waves <- c("one","two","three","four","five")
locationsYP <- sapply(waves, function(w) paste0("./UKDA-5545-spss/spss/spss28/eul/wave_",w,"_lsype_young_person_2020.sav")  )
locationsYPList <- as.list(locationsYP)
dataListYP <- lapply(locationsYPList, read.spss, use.value.labels = TRUE, to.data.frame = TRUE)

#remove approx 1000 observations where YP surveyed in September and October in Wave 5 since they will have finished Year 13 and not be eligible for EMA 
dataListYP$five <- dataListYP$five %>% filter(!(W5IntMonth %in% c("September","October")))

locationsHH <- sapply(waves, function(w) paste0("./UKDA-5545-spss/spss/spss28/eul/wave_",w,"_lsype_family_background_2020.sav")  )
locationsHHList <- as.list(locationsHH)
dataListHH <- lapply(locationsHHList, read.spss, use.value.labels = TRUE, to.data.frame = TRUE)

dataJoin <- dataListYP$three %>% 
  left_join(dataListHH$one, by=join_by(NSID==NSID)) %>% 
  left_join(dataListHH$two, by=join_by(NSID==NSID)) %>% 
  left_join(dataListHH$three, by=join_by(NSID==NSID)) %>% 
  left_join(dataListYP$one, by=join_by(NSID==NSID)) %>%
  left_join(dataListYP$two, by=join_by(NSID==NSID)) %>%
  left_join(dataListYP$four, by=join_by(NSID==NSID)) %>% 
  left_join(dataListYP$five, by=join_by(NSID==NSID)) %>%
  rename(W2gor=gor.x, W3gor=gor.y, W1ethgrpYP=W1ethgrpYP.x)

W4nsid <- dataListYP$four$NSID
W5nsid <- dataListYP$five$NSID

data <- dataJoin %>% 
  mutate(W1GrssyrHH=as.numeric(as.character(W1GrssyrHH)),
         W2GrssyrHH=as.numeric(as.character(W2GrssyrHH)),
         W4EMA2YP=as.numeric(as.character(W4EMA2YP)),
         W5EMA2YP=as.numeric(as.character(W5EMA2YP)),
         W2BenSalAnTot=as.numeric(as.character(W2BenSalAnTot)),
         W3BenSalAnTotUpper=as.numeric(gsub("[£,]", "",str_extract(W3incestm, "£[0-9,]+$")))*12,
         W3BenSalAnTotLower=as.numeric(gsub("[£,]", "",str_extract(W3incestm, "£[0-9,]+")))*12) %>% 
  mutate(W3BenSalAnTotLower=ifelse(W3BenSalAnTotLower==2592,0,W3BenSalAnTotLower),
         W3BenSalAnTotUpper=ifelse(W3BenSalAnTotLower==51996,100000,W3BenSalAnTotUpper)) %>%
  mutate(W4EMA2YP=ifelse(is.na(W4EMA2YP),0,W4EMA2YP),
         W5EMA2YP=ifelse(is.na(W5EMA2YP),0,W5EMA2YP)) %>%
  mutate(W4cens = !(NSID %in% W4nsid),
         W5cens = (NSID %in% W4nsid) & !(NSID %in% W5nsid),
         W4edu = (W4MainActYP == "Going to a school or college full time "),
         W5edu = (W5mainactYP == "In education"),
         W4highEMA = ifelse(!is.na(W4EMA2YP) &  W4EMA2YP >= 30,1,0)) %>%
  filter(!(W4cens==FALSE & is.na(W4edu))) %>%
  filter(!(W4cens==FALSE & W5cens==FALSE & is.na(W5edu)))

set.seed(1)
data$u <- runif(n=nrow(data))
data <- data %>% 
  mutate(W3SalEst = case_when(
    (!is.na(is.na(W3incestm)) & W2BenSalAnTot <= W3BenSalAnTotUpper & W2BenSalAnTot >= W3BenSalAnTotLower) | (is.na(W3incestm) & !is.na(W2BenSalAnTot)) ~ W2BenSalAnTot,
    (!is.na(is.na(W3incestm)) & W3BenSalAnTotLower>=51996 & !is.na(W2BenSalAnTot) & W2BenSalAnTot >= 50000) ~ W2BenSalAnTot,
    (!is.na(is.na(W3incestm)) & W3BenSalAnTotLower>=51996 & !is.na(W2GrssyrHH) & W2GrssyrHH >= 50000) ~ W2GrssyrHH,
    (!is.na(is.na(W3incestm)) & W3BenSalAnTotLower>=51996 & !is.na(W1GrssyrHH) & W1GrssyrHH >= 50000) ~ W1GrssyrHH,
    .default = u*W3BenSalAnTotUpper+(1-u)*W3BenSalAnTotLower #Note: The W3incestm intervals luckily line up almost exactly with the EMA thresholds, so we don't risk moving people across the cutoffs
  )) %>% 
  filter(!is.na(W3SalEst))

D1ids <- sample(1:nrow(data), size=round(nrow(data)/2), replace=FALSE)
data[D1ids,"D1"] <- 1
data[-D1ids,"D1"] <- 0

EMAIncomeThresholds <- c(20817,25521,30810)

##Step 0: Initialize data frames for estimation
outW4 <- c("NSID","D1","W4cens","W4edu","W4highEMA","W5cens","W5edu")
outW5 <- outW4[!(outW4 %in% c("W4cens","W4highEMA"))]

dfFull <- data %>% 
  select(any_of(outW4),starts_with("W1"),starts_with("W2"),starts_with("W3")) 

df <- data %>% 
  select(any_of(outW4),W3SalEst,
         W3hous12HH,W3usevcHH,W3ageMP,W3famtyp,
         W3sexMP,W3empsMP,W3cnsseccatfam,W3gor,
         W3depkids,W3sexYP,W3jobYP,W3jobearnYP,W1ethgrpYP) 

D1 <- df %>% filter(D1==1)
D2 <- df %>% filter(D1==0)

rmOutW4 <- function(df){
  df %>% select(-any_of(outW4))
}
rmOutW5 <- function(df){
  df %>% select(-any_of(outW5))
}
rmW4Cens <- function(df){
  df %>% filter(W4cens==FALSE)
}
rmW5Cens <- function(df){
  df %>% filter(W4cens==FALSE & W5cens==FALSE)
}


##Step 1.1: Estimate censoring distribution

#Find hyperparameters
paramsCens <- list(
  objective = "binary:logistic",
  eval_metric = "logloss", 
  eval_metric = "auc",
  eta = 0.05,
  max_depth = 5,
  min_child_weight = 5,
  subsample = 0.6,
  colsample_bytree = 0.7)

paramsCensFull <- list(
  objective = "binary:logistic", 
  eval_metric = "logloss",
  eval_metric = "auc",
  eta = 0.01,
  max_depth = 20, 
  min_child_weight = 10,
  subsample = 0.5,
  colsample_bytree = 0.5,
  scale_pos_weight = 0.8*sum(data$W4cens==0)/sum(data$W4cens==1),
  alpha=0.1)

set.seed(1)
#xgb.cv(params=paramsCens,data=data.matrix(rmOutW4(df)), label=dfR$W4cens, nfold=5, nrounds=1000, early_stopping_rounds = 200)
#xgb.cv(params=paramsCensFull,data=data.matrix(rmOutW4(dfFull)), label=dfFull$W4cens, nfold=5, nrounds=1000, early_stopping_rounds = 200) 
gc()

set.seed(1)
#xgb.cv(params=paramsCens,data=data.matrix(rmOutW5(rmW4Cens(df))), label=rmW4Cens(df)$W5cens, nfold=5, nrounds=1000, early_stopping_rounds = 200) 
#xgb.cv(params=paramsCensFull,data=data.matrix(rmOutW5(rmW4Cens(dfFull))), label=rmW4Cens(dfFull)$W5cens, nfold=5, nrounds=1000, early_stopping_rounds = 200)
gc()


#Wave 4
set.seed(1)
CensD1W4 <- xgboost(data=data.matrix(rmOutW4(D1)), params=paramsCens, label=D1$W4cens,nrounds=1000, early_stopping_rounds = 200)
CensD2W4 <- xgboost(data=data.matrix(rmOutW4(D2)), params=paramsCens, label=D2$W4cens,nrounds=1000, early_stopping_rounds = 200)

#Wave 5
set.seed(1)
CensD1W5 <- xgboost(data=data.matrix(rmOutW5(rmW4Cens(D1))), params=paramsCens, label=rmW4Cens(D1)$W5cens, nrounds=1000, early_stopping_rounds = 200)
CensD2W5 <- xgboost(data=data.matrix(rmOutW5(rmW4Cens(D2))), params=paramsCens, label=rmW4Cens(D2)$W5cens, nrounds=1000, early_stopping_rounds = 200)

##Step 1.2: Estimate outcome distribution

#find hyperparameters
paramsOut <- list(
  objective = "binary:logistic",
  eval_metric = "logloss", 
  eval_metric = "auc",
  eta = 0.1,
  max_depth = 5,
  min_child_weight = 5, 
  subsample = 0.7,
  colsample_bytree = 0.7)

paramsOutFull <- list(
  objective = "binary:logistic",
  eval_metric = "logloss", 
  eval_metric = "auc",
  eta = 0.06,
  max_depth = 10,
  min_child_weight = 20,
  subsample = 0.4,
  colsample_bytree = 0.5)

set.seed(1)
#xgb.cv(params=paramsOut,data=data.matrix(rmOutW4(rmW4Cens(df))), label=rmW4Cens(df)$W4edu, nfold=5, nrounds=1000, early_stopping_rounds = 200)
#xgb.cv(params=paramsOutFull,data=data.matrix(rmOutW4(rmW4Cens(dfFull))), label=rmW4Cens(dfFull)$W4edu, nfold=5, nrounds=1000, early_stopping_rounds = 200)
gc()

set.seed(1)
#xgb.cv(params=paramsOut,data=data.matrix(rmOutW5(rmW5Cens(df))), label=rmW5Cens(df)$W5edu, nfold=5, nrounds=1000, early_stopping_rounds = 200)
#xgb.cv(params=paramsOutFull,data=data.matrix(rmOutW5(rmW5Cens(dfFull))), label=rmW5Cens(dfFull)$W5edu, nfold=5, nrounds=1000, early_stopping_rounds = 200)
gc()

#Wave 4
set.seed(1)
OutcomeD1W4 <- xgboost(data=data.matrix(rmOutW4(rmW4Cens(D1))), params=paramsOut, label=rmW4Cens(D1)$W4edu,nrounds=1000, early_stopping_rounds = 200)
OutcomeD2W4 <- xgboost(data=data.matrix(rmOutW4(rmW4Cens(D2))), params=paramsOut, label=rmW4Cens(D2)$W4edu,nrounds=1000, early_stopping_rounds = 200)

#Wave 5
set.seed(1)
OutcomeD1W5 <- xgboost(data=data.matrix(rmOutW5(rmW5Cens(D1))), params=paramsOut, label=rmW5Cens(D1)$W5edu,nrounds=1000, early_stopping_rounds = 200)
OutcomeD2W5 <- xgboost(data=data.matrix(rmOutW5(rmW5Cens(D2))), params=paramsOut, label=rmW5Cens(D2)$W5edu,nrounds=1000, early_stopping_rounds = 200)

##Step 2.1: Create pseudo-outcomes

pseudoValues <- function(df){
  
  if(df$D1[1]){
    CensW5 <- CensD2W5
    CensW4 <- CensD2W4
    OutcomeW5 <- OutcomeD2W5
    OutcomeW4 <- OutcomeD2W4
  }
  else{
    CensW5 <- CensD1W5
    CensW4 <- CensD1W4
    OutcomeW5 <- OutcomeD1W5
    OutcomeW4 <- OutcomeD1W4
  }
  
  dfPseudo <- df
  
  #Compute censoring probabilities
  dfPseudotmp <- df %>% rmW4Cens() 
  dfPseudotmp$pCW5 <- predict(object=CensW5, newdata=data.matrix(rmOutW5(dfPseudotmp)))
  dfPseudotmp <- dfPseudotmp %>% select(NSID,pCW5)
  
  dfPseudo$pCW4 <- predict(object=CensW4, newdata=data.matrix(rmOutW4(df)))
  dfPseudo <- dfPseudo %>% left_join(dfPseudotmp, by=join_by(NSID==NSID))
  
  #Compute conditional expectations
  dfPseudotmp0 <- df %>% mutate(W4edu=0)
  dfPseudotmp0$EduW5NoEduW4 <- predict(object=OutcomeW5, newdata=data.matrix(rmOutW5(dfPseudotmp0)))  
  dfPseudotmp0 <- dfPseudotmp0 %>% select(NSID,EduW5NoEduW4)
  
  dfPseudotmp1 <- df %>% mutate(W4edu=1)
  dfPseudotmp1$EduW5EduW4 <- predict(object=OutcomeW5, newdata=data.matrix(rmOutW5(dfPseudotmp1)))     
  dfPseudotmp1 <- dfPseudotmp1 %>% select(NSID,EduW5EduW4)
  
  dfPseudotmp <- df 
  dfPseudotmp$EduW4 <- predict(object=OutcomeW4, newdata=data.matrix(rmOutW4(dfPseudotmp)))     
  dfPseudotmp <- dfPseudotmp %>% select(NSID,EduW4)
  
  dfPseudo <- dfPseudo %>% 
    left_join(dfPseudotmp0, by=join_by(NSID==NSID)) %>% 
    left_join(dfPseudotmp1, by=join_by(NSID==NSID)) %>%
    left_join(dfPseudotmp, by=join_by(NSID==NSID))
  
  dfPseudo <- dfPseudo %>% 
    mutate(VY1 = EduW4*1 + (1-EduW4)*0 + EduW4*EduW5EduW4 + (1-EduW4)*EduW5NoEduW4,
           VY2 = W4edu*EduW5EduW4 + (1-W4edu)*EduW5NoEduW4 )
  
  dfPseudo <- dfPseudo %>%
    mutate(integral1 = VY1/(1-pCW4)*(W4cens-pCW4),
           integral2 = (W4edu+VY2)/((1-pCW4)*(1-pCW5))*(W5cens-(1-W4cens)*pCW5),
           IPCW = ifelse(W4cens==0 & W5cens==0,W4edu+W5edu,0)/((1-pCW4)*(1-pCW5)),
           IPCWA = ifelse(W4cens==0,W4highEMA/(1-pCW4),0) ) %>%
    mutate(integral1=ifelse(is.na(integral1),0,integral1),
           integral2=ifelse(is.na(integral2),0,integral2),
           IPCW=ifelse(is.na(IPCW),0,IPCW),
           IPCWA=ifelse(is.na(IPCWA),0,IPCWA))
  
  dfPseudo <- dfPseudo %>% mutate(Ystar = IPCW+integral1+integral2)
  
  dfPseudo
}

D1Pseudo <- pseudoValues(D1)
D2Pseudo <- pseudoValues(D2)

D1Pseudo <- D1Pseudo %>% filter(!is.na(W3SalEst))
D2Pseudo <- D2Pseudo %>% filter(!is.na(W3SalEst))

##Step 2.2: Pseudo-outcome regression

b <- 3500 #4704 #2352 

fit1 <- rdrobust(x=D1Pseudo$W3SalEst,y=D1Pseudo$Ystar,fuzzy=D1Pseudo$IPCWA, c=EMAIncomeThresholds[1], h=b)
fit2 <- rdrobust(x=D2Pseudo$W3SalEst,y=D2Pseudo$Ystar,fuzzy=D2Pseudo$IPCWA, c=EMAIncomeThresholds[1], h=b)

fit1$N_h
fit2$N_h

##Step 3: Cross-fitting

SeConstant <- 1/sqrt(2)
fitCross <- fit1$Estimate[,"tau.us"]*0.5+fit2$Estimate[,"tau.us"]*0.5
SECross <-  SeConstant*(fit1$Estimate[,"se.us"]*0.5+fit2$Estimate[,"se.us"]*0.5)

fitCross 
SECross

##Plots
DfPseudoPlot <- union(D1Pseudo,D2Pseudo) %>% 
  mutate(EMAtreat = (W3SalEst <= EMAIncomeThresholds[1])) 

#Outcome global
DfPseudoPlotNoRD <- DfPseudoPlot
DfPseudoPlotNoRD$bins <- cut(DfPseudoPlotNoRD$W3SalEst, breaks = 400)
binned_meansNoRD <- aggregate(Ystar ~ bins, data = DfPseudoPlotNoRD, FUN = mean)
binned_meansNoRD$midpoint <- DfPseudoPlotNoRD %>% group_by(bins) %>% summarize(W3SalEst=mean(W3SalEst)) %>% select(W3SalEst) %>% unlist()
pSal <- ggplot(binned_meansNoRD, aes(x = midpoint, y = Ystar), color="black") +
  geom_step() +
  labs(x = "Salary", y = "Mean of pseudo-outcome") +
  theme_light() + 
  geom_vline(xintercept = EMAIncomeThresholds[1], linetype="dotted", 
             color = "gray50", size=1.5) +
  geom_vline(xintercept = EMAIncomeThresholds[2], linetype="dotted", 
             color = "gray50", size=1.5) +
  geom_vline(xintercept = EMAIncomeThresholds[3], linetype="dotted", 
             color = "gray50", size=1.5) +
  xlim(c(0,1.1*10^5)) + 
  ylim(c(1.1,2)) +
  theme(plot.title = element_text(size=18),
        text = element_text(size = 18),
        strip.text = element_text(size=18, face="bold"),
        legend.text=element_text(size=16)) 

#Outcome close to threshold
grpK <- function(df,k){
  df %>% filter(EMAtreat == k) %>% filter(W3SalEst <= EMAIncomeThresholds[2])
}

reg0 <- lprobust(y=grpK(DfPseudoPlot,0)$Ystar,x=grpK(DfPseudoPlot,0)$W3SalEst,eval=seq(from=EMAIncomeThresholds[1],to=EMAIncomeThresholds[2],by=200),h=b)$Estimate[,c("eval","tau.us","se.us")] %>% 
  as.data.frame() %>% 
  mutate(CIlower = tau.us-1.96*se.us, CIupper = tau.us+1.96*se.us)
reg1 <- lprobust(grpK(DfPseudoPlot,1)$Ystar,grpK(DfPseudoPlot,1)$W3SalEst,eval=seq(from=0,to=EMAIncomeThresholds[1],by=200),h=b)$Estimate[,c("eval","tau.us","se.us")] %>% 
  as.data.frame() %>%
  mutate(CIlower = tau.us-1.96*se.us, CIupper = tau.us+1.96*se.us)

DfPseudoPlotRD <- DfPseudoPlot %>% filter(W3SalEst <= EMAIncomeThresholds[2])
DfPseudoPlotRD$bins <- cut(DfPseudoPlotRD$W3SalEst, breaks = 100) #100
binned_meansRD <- aggregate(Ystar ~ bins, data = DfPseudoPlotRD, FUN = mean)
binned_meansRD$midpoint <- DfPseudoPlotRD %>% group_by(bins) %>% summarize(W3SalEst=mean(W3SalEst)) %>% select(W3SalEst) %>% unlist()
pOut <- ggplot(binned_meansRD, aes(x = midpoint, y = Ystar, color="black")) +
  geom_step() +
  theme_light() +
  labs(x = "Salary", y = "Mean of pseudo-outcome") +
  geom_vline(xintercept = EMAIncomeThresholds[1], linetype="dotted", 
             color = "gray50", size=2) +
  geom_vline(xintercept = EMAIncomeThresholds[2], linetype="dotted", 
             color = "gray50", size=2) +
  geom_vline(xintercept = EMAIncomeThresholds[3], linetype="dotted", 
             color = "gray50", size=2) +
  xlim(c(15000,EMAIncomeThresholds[2])) + 
  ylim(c(1,1.65)) +
  theme(plot.title = element_text(size=18),
        text = element_text(size = 18),
        strip.text = element_text(size=18, face="bold"),
        legend.text=element_text(size=16),
        legend.position = c(0.29,0.92),
        legend.background = element_rect(fill='transparent'),
        legend.key = element_rect(fill = "transparent"))
pOut <- pOut + 
  geom_line(data=reg0, aes(x=eval, y = tau.us, color="ForestGreen"),lwd=1.1) +
  geom_line(data=reg1, aes(x=eval, y = tau.us, color="ForestGreen"),lwd=1.1) +
  geom_line(data=reg0, aes(x=eval, y = CIupper, color="blue"), lty=2, lwd=1) + 
  geom_line(data=reg0, aes(x=eval, y = CIlower, color="blue"), lty=2, lwd=1) + 
  geom_line(data=reg1, aes(x=eval, y = CIupper, color="blue"), lty=2, lwd=1) + 
  geom_line(data=reg1, aes(x=eval, y = CIlower, color="blue"), lty=2, lwd=1) +
  scale_colour_manual(name = "", labels = c("Binned means", "95% CI", "Estimate"), values = c("black", "blue", "ForestGreen"))
  
#Treatment
regA0 <- lprobust(grpK(DfPseudoPlot,0)$W4highEMA,grpK(DfPseudoPlot,0)$W3SalEst,eval=seq(from=EMAIncomeThresholds[1],to=EMAIncomeThresholds[2],by=200),h=b)$Estimate[,c("eval","tau.us","se.us")] %>% 
  as.data.frame() %>% 
  mutate(CIlower = tau.us-1.96*se.us, CIupper = tau.us+1.96*se.us)
regA1 <- lprobust(grpK(DfPseudoPlot,1)$W4highEMA,grpK(DfPseudoPlot,1)$W3SalEst,eval=seq(from=0,to=EMAIncomeThresholds[1],by=200),h=b)$Estimate[,c("eval","tau.us","se.us")] %>% 
  as.data.frame() %>% 
  mutate(CIlower = tau.us-1.96*se.us, CIupper = tau.us+1.96*se.us)

DfPseudoPlot$bins <- cut(DfPseudoPlot$W3SalEst, breaks = 1600)
binned_meansA <- aggregate(W4highEMA ~ bins, data = DfPseudoPlot, FUN = mean)
binned_meansA$midpoint <- DfPseudoPlot %>% group_by(bins) %>% summarize(W3SalEst=mean(W3SalEst)) %>% select(W3SalEst) %>% unlist()
binned_meansA <- binned_meansA
pA <- ggplot(binned_meansA, aes(x = midpoint, y = W4highEMA, color="black")) +
  geom_step() +
  labs(x = "Salary", y = "Probability of receiving EMA") +
  theme_light() + 
  geom_vline(xintercept = EMAIncomeThresholds[1], linetype="dotted", 
             color = "gray50", size=2) +
  geom_vline(xintercept = EMAIncomeThresholds[2], linetype="dotted", 
             color = "gray50", size=2) +
  geom_vline(xintercept = EMAIncomeThresholds[3], linetype="dotted", 
             color = "gray50", size=2) +
  xlim(c(15000,EMAIncomeThresholds[2])) +
  ylim(c(0.2,0.7)) +
  theme(plot.title = element_text(size=18),
        text = element_text(size = 18),
        strip.text = element_text(size=18, face="bold"),
        legend.text=element_text(size=16),
        legend.position = c(0.29,0.92),
        legend.background = element_rect(fill='transparent'),
        legend.key = element_rect(fill = "transparent"))
pA <- pA + 
  geom_line(data=regA0, aes(x=eval, y=tau.us, color="ForestGreen"), size=1.1) +
  geom_line(data=regA1, aes(x=eval, y=tau.us, color="ForestGreen"), size=1.1) +
  geom_line(data=regA0, aes(x=eval, y=CIlower, color="blue"),lwd=1.1, lty=2) +
  geom_line(data=regA0, aes(x=eval, y=CIupper, color="blue"),lwd=1.1, lty=2) +
  geom_line(data=regA1, aes(x=eval, y=CIlower, color="blue"),lwd=1.1, lty=2) +
  geom_line(data=regA1, aes(x=eval, y=CIupper, color="blue"),lwd=1.1, lty=2) +
  scale_colour_manual(name = "", labels = c("Binned means", "95% CI", "Estimate"), values = c("black", "blue", "ForestGreen"))

pAll <- grid.arrange(pSal, pOut, pA, nrow = 1)

ggsave(filename="Figures/RDDPlot.png",plot=pAll,height=6, width=14)
