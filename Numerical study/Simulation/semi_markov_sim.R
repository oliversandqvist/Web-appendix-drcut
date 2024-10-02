#Reference: "On Lewis' simulation method for point processes". For more details, "Exact simulation of the Jump Times of a Class of Piecewise Deterministic Markov Processes".

no_states <- 4
Xmax <- 4
nmax <- 3

#Computes the transition rate for a given transition at a given time and duration
#x is current state #y is destination state #s is time of last jump #t is current time #cov is the baseline covariates known at time 0 (we use age0 = age at time 0 and gender is gender where female is 0 and male is 1)
mu12 <- function(t,cov){
  return(exp(log(0.3)+0.15*cos(0.5*pi*cov)+0.15*(t > 2.5)-0.05*cov)) 
}

mu13 <- function(t,cov){
  return(exp(log(0.1)+0.3*sin(0.5*pi*cov)+0.05*t))
}

mu23 <- function(t,s,cov){
  covPol <- min(cov,3)
  return(exp(-0.75*min(t-s,3)*(1.07+0.09*covPol-0.024*covPol^2-0.014*covPol^3+0.001*covPol^4+0.00065*covPol^5)))
}

mu14 <- function(t,s,cov){
  return(exp(log(0.2)+0.6*(cov >= -2 & cov < 2))) 
}

mu24 <- function(t,s,cov){
  return(0) 
}

transition_rate <-  function(x, y, s, t, cov){
  if(x == 1 && y==2){
    as.numeric(mu12(t,cov))
  }
  else if(x == 1 && y==3){
    as.numeric(mu13(t,cov)) 
  }
  else if(x==1 && y==4){
    as.numeric(mu14(t,s,cov))
  }
  else if(x==2 && y==3){
    as.numeric(mu23(t,s,cov))
  }
  else if(x==2 && y==4){
    as.numeric(mu24(t,s,cov))
  }
  else{0}
}


#Computes the total intensity of out the current state
#x is current state #s is time of last jump #t is current time #cov is the baseline covariates known at time 0
jump_rate <-  function(x, s, t, cov){
  sum(sapply(1:no_states, FUN = function(y) transition_rate(x,y,s,t,cov)))
}

#Computes the mark distribution given a jump time
#x is current state #s is time of last jump #t is time of current jump
mark_dist <- function(x, s, t, cov){
  rates <- sapply(1:no_states, FUN = function(y) transition_rate(x,y,s,t,cov))
  rates/sum(rates)
}

#Simulates time from the last jump to the next jump 
#x is current state #s is time of last jump #b is bound #cov is the baseline covariates known at time 0
jump_simulator <- function(x, s, tstop, b, cov){ 
  u <- runif(1)
  e <- rexp(1, rate = 1)
  t <- e/b
  if (t > tstop) { return(Inf) }
  while(u > jump_rate(x, s, s+t, cov)/b){ 
    u <- runif(1)
    e <- rexp(1, rate = 1)
    t <- e/b+t
    if (t > tstop) { return(Inf) }
  }
  return(t)
}

#Computes the "local upper bounds" of the intensities for each of the states.
#s is the time of the last jump at time 0 #t is the current time #tstop is the terminal time #cov is the baseline covariates known at time 0
rate_bounds <- function(s, t, tstop, cov){
  b12 <- transition_rate(1,2,0,tstop,cov)
  b13 <- transition_rate(1,3,0,tstop,cov)
  b14 <- transition_rate(1,4,0,tstop,cov)
  b1 <- b12+b13+b14
  b23 <- transition_rate(2,3,0,0,cov)
  b24 <- transition_rate(2,4,tstop,tstop,cov)
  b2 <- b23+b24
  b3 <- 0
  b4 <- 0
  rate_bounds <- c(b1,b2,b3,b4)
  return(rate_bounds)
}

#Simulates one trajectory of the semi-Markov process.
#x0 is initial state #s0 is last jump time #t0 is current time (note that by convention, we include the zero'th jump-time and mark in the event history) #tstop is terminal time of the trajectory #cov is the baseline covariates known at time 0
simulator <- function(x0, s0, t0, tstop, cov,AJformat=FALSE){ 
  rate_bounds <- rate_bounds(s0, t0, tstop, cov)
  
  times <- c(t0)
  states <- c(x0)
  lastJump <- s0
  lastState <- x0
  t <- t0
  bx <- rate_bounds[lastState]
  
  if(bx == 0) { break } #break if absorbed
  
  repeat{
    t <- lastJump+jump_simulator(lastState, lastJump, tstop, bx, cov)
    
    if(t > tstop){ break } #break if jump time exceeds terminal time
    
    x <- sample(1:no_states, 1, prob = mark_dist(lastState, lastJump, t, cov))
    
    times <- c(times, t)
    states <- c(states, x)
    lastJump <- t
    lastState <- x
    bx <- rate_bounds[lastState]
    
    if(bx == 0){ break } #break if absorbed
    
  }
  
  if(AJformat){
    if(bx != 0){ #for transient states, the last state and jump is repeated at the censoring time
      times <- c(times, tstop)
      states <- c(states, lastState)
    }
    if(lastState==4){ #censorings are not to be included as transitions as this affects transition probabilities
      len <- length(states)
      states[len] <- states[len-1] 
    }
  }
  
  return(list(times = times, states = states))
}

simulator_n <- function(x0, s0, t0, tstop, cov, AJformat=FALSE){ 
  n <- length(cov)
  namesData <- c("id","X","tn","tnNext","yn","ynNext")
  nCol <- length(namesData)
  simDf <- setNames(data.table(matrix(nrow = n*nmax, ncol = nCol)), namesData)
  simDf <- simDf[, id:=as.integer(id)]
  simDf <- simDf[, lapply(.SD, as.numeric), by=id]
  simList <- vector("list", n) #list()
  
  k <- 1L
  for(m in 1:n){
    covm <- cov[m]
    rate_bounds <- rate_bounds(s0, t0, tstop, covm)
    
    times <- c(t0)
    states <- c(x0)
    lastJump <- s0
    lastState <- x0
    t <- t0
    bx <- rate_bounds[lastState]
    
    if(bx == 0) { break } #break if absorbed
    
    repeat{
      t <- lastJump+jump_simulator(lastState, lastJump, tstop, bx, covm)
      
      if(t > tstop){ break } #break if jump time exceeds terminal time
      
      x <- sample(1:no_states, 1, prob = mark_dist(lastState, lastJump, t, covm))
      
      times <- c(times, t)
      states <- c(states, x)
      
      if(!AJformat){
        set(simDf, i=k, j=1L, value=m)
        set(simDf, i=k, j=2L, value=covm)
        set(simDf, i=k, j=3L, value=lastJump)
        set(simDf, i=k, j=4L, value=t)    
        set(simDf, i=k, j=5L, value=lastState)   
        set(simDf, i=k, j=6L, value=x)  
        k <- k+1L
      }
      
      lastJump <- t
      lastState <- x
      bx <- rate_bounds[lastState]
      
      
      if(bx == 0){ break } #break if absorbed
      
    }
    
    if(AJformat){
      if(bx != 0){ #for transient states, the last state and jump is repeated at the censoring time
        times <- c(times, tstop)
        states <- c(states, lastState)
      }
      if(lastState==4){ #censorings are not to be included as transitions as this affects transition probabilities
        len <- length(states)
        states[len] <- states[len-1] 
      }
      simList[[m]]$times <- times
      simList[[m]]$states <- states
      simList[[m]]$X <- covm
    }
    else{
      set(simDf, i=k, j=1L, value=m)
      set(simDf, i=k, j=2L, value=covm)
      set(simDf, i=k, j=3L, value=lastJump)
      set(simDf, i=k, j=4L, value=NA)    
      set(simDf, i=k, j=5L, value=lastState)   
      set(simDf, i=k, j=6L, value=NA)  

      k <- k+1L
    }
  }
  
  if(AJformat){
    return(simList)
  }
  else{
    return(simDf)
  }
}
