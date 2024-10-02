##' @param covars covariates to include in HAL model. 
##' @param dt dataset.
##' @param V number of folds in cross-validation.
##' @param cut number of cuts for each covariate.
##' @param seed seed used for cross-validation.
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here-
##' @param time.var name of time variable (note: if time.var is included in fixed.covars, it will still increase with time)
fit.hal <- function(fixed.covars, time.covars, dt, V=5, cut=8, seed=13349, 
                    delta.var="delta", delta.value=1,
                    time.var=NULL,browse1=FALSE) {
    
    set.seed(seed)
  
    if (browse1) browser()
    
    if (length(time.var)==0){error("error: no time variable.")}
  
    X.hal <- basis.fun(fixed.covars=fixed.covars, time.covars=time.covars, 
                            cut=cut, dt=dt,
                            time.var=time.var,
                            delta.var=delta.var, delta.value=delta.value)
    
    pseudo.dt <- X.hal$pseudo.dt
    model.formula.X <- X.hal$model.formula.X 
    X.hal <- X.hal$X
    
    pseudo.dt[, D:=sum(time.obs==grid.time.right & get(delta.var)==delta.value), by="x"] #find occurrences for the transition of interest (delta.value) for a given combination of indicator functions in HAL (x)
    pseudo.dt[, RT:=sum(risk.time), by="x"] #find exposures for a given combination of indicator functions in HAL (x); note that the code is in a "competing risk" form, so don't have to distinguish between which state the exposure pertains to. For intensity estimation, therefore have to partition data into a competing risk form before feeding into fit.hal 
    
    tmp.dt <- unique(pseudo.dt[, c("x", "D", "RT")]) #each id with a given combination of indicator functions from HAL (x) will have a copy of the corresponding value of occurrences and exposures for this group; here we only take the unique values (for each x)
    Y2.hal <- tmp.dt[RT>0, D]
    offset2 <- tmp.dt[RT>0, log(RT)]
    X2.hal <- unique.matrix(X.hal)[tmp.dt$RT>0,]
    
    #fit HAL with glmnet
    cve.hal <- cv.glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                      offset=offset2,
                      family=poisson(link="log"),
                      maxit=100000,
                      nfolds=V)
    hal.fit <- cve.hal$glmnet.fit
    
    return(list(fit=hal.fit,
                cve=cve.hal,
                model.formula.X=model.formula.X))
}
