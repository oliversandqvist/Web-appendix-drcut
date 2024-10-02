##' @param covars covariates to include in HAL model. 
##' @param cut number of cuts for each covariate.
##' @param dt dataset. 
##' @param Y outcome object. 
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here.
##' @param time.var name of time-variable. 
basis.fun <- function(fixed.covars, time.covars, cut, dt,
                      Y=paste0("(", delta.var, "==", delta.value, ")"),
                      delta.var="delta", delta.value=1,
                      time.var=NULL) {
  
    covars <- c(fixed.covars,time.covars)
    if (length(covars)>0) { #make covariate columns to numeric type
      dt[, (covars):=lapply(.SD, function(x) {
        if (is.character(x)) {
          return(as.numeric(as.factor(x)))
        } else return(x)
      }), .SDcols=covars]
    }
  
    pseudo.dt <- copy(dt)        
    grid.times <- c(0, indicator.basis.fun(dt, time.var, cut, return.grid=TRUE))
    pseudo.dt <- do.call("rbind", lapply(1:length(grid.times), function(jj) { #copy each row of pseudo.dt and attaches each of the time cut points to each of these new rows
        tmp <- copy(pseudo.dt)[, grid.period:=jj] 
    }))
    pseudo.dt[, grid.time.left:=grid.times[grid.period]] #here, grid.time is left end-point of interval
    pseudo.dt <- pseudo.dt[grid.time.left<get(time.var)] #filter away the rows where the left end-point is larger than the event time 
    pseudo.dt <- pseudo.dt[order(id, get(time.var))]
    
    #update time.covars linearly with passing of time
    if (length(time.covars)>0) { 
      for(cov in time.covars){
        pseudo.dt <- pseudo.dt[, (cov):=get(cov)+grid.time.left]
      }
    }
    
    pseudo.dt[, grid.time.right:=grid.times[grid.period+1]] #increases time.grid to the next cut point so it is right end-point of interval
    pseudo.dt[get(time.var)<=grid.time.right, grid.time.right:=get(time.var)] #cap grid.time at time observation (so we get grid.time intersected [0,time])
    pseudo.dt[, time.obs:=get(time.var)] #saves the time observations in time.obs column, and in the following line, we overwrite time observations by time.grid
    pseudo.dt[, (time.var):=grid.time.left]
    

    indicator.basis.list <- lapply(covars,FUN=function(cov) indicator.basis.fun(pseudo.dt,cov,cut) )
    model.formula.X <- formula(paste0(
      "~-1",
      paste0("+",paste0(apply(do.call(what=expand.grid,args = indicator.basis.list), 1, function(row) paste0(row, collapse="*")),collapse="+"))
    ))
    model.formula <- formula(paste0(Y, paste0(model.formula.X, collapse=""), collapse = ""))
    X <- Matrix(model.matrix(model.formula, data=pseudo.dt), sparse=FALSE)
    
    x.vector <- apply(X, 1, function(x) paste0(x, collapse=",")) #collapse matrix to vector where each element of vector is the string of 0's and 1's of the columns of the matrix for that row #NB: This is the slowest part of the code
    pseudo.dt[, x:=x.vector]
    pseudo.dt[, risk.time:=grid.time.right-grid.time.left, by="id"] #for each id, take the difference in grid.time to get the exposure in that interval
    return(list(X=X, pseudo.dt=pseudo.dt, model.formula.X = model.formula.X)) #return design matrix and dataframe modified with times used for exposures
}

indicator.basis.fun <- function(mat, xvar, xcut, type="obs", return.grid=FALSE, seed=13349) { #cut sort values of covariate "xvar" and cut into "xcut" pieces. Return a vector of indicator events, where the n'th event is that xvar exceeds the n'th cut-point
    set.seed(seed)
    if (type=="obs") {
        xvar.values <- mat[, sort(unique(get(xvar)))]
        xvar.pick <- seq(1, length(xvar.values), length=min(xcut, length(xvar.values)))
        xgrid <- xvar.values[xvar.pick][-c(1)]
    } else {
        xgrid <- round(seq(mat[, min(get(xvar))],
                           mat[, max(get(xvar))],
                           length=xcut)[-c(1)], 2)
    }
    if (return.grid) return(c(unique(xgrid))) else return(paste0("(", xvar, ">=", unique(xgrid), ")"))
}

