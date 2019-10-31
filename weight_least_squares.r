lm.wls <- function(y, X, V=rep(1, length(y))) {
  ####################################################################################
  ##
  ##    A function that calculates the weighted least squares (WLS) estimates for a 
  ##    regression model given a dataset and a vector of weights. 
  ##
  ##    The arguments are:
  ##
  ##    y : An n-vector of the response variables of the model.
  ##    X : An nxp matrix of the covariate values for each respective y-value.
  ##    V : An n-vector of weights for each covariate in the model. Defaults to a
  ##        vector of 1s, in which case regular least squares estimation occurs.
  ##
  ##    The function returns a named list containing:
  ##    
  ##    betahat : A p-vector of the WLS estimates for the model covariates. 
  ##    sebeta  : A p-vector of the standard errors for the WLS estimates. 
  ##    sigsq   : The estimated variance for each observation.
  ##
  ####################################################################################
  # 
  #  Start by defining the model's dimensions n and p, the number of observations
  #  and covariates respectively.
  #
  n <- length(y); p <- ncol(X)  
  #
  #  Calculate betahat, the WLS estimator of the model coefficients. We use a tryCatch
  #  block to catch errors relating to the arguments' dimensions and provide our own
  #  more helpful error messages. We do this inside a withCallingHandlers block
  #  to suppress unnecessary warning messages arising after our errors have been thrown.
  #
  V <- 1/V       # Redefine V as a vector of its reciprocals for use later on.
  betahat <- withCallingHandlers(
    tryCatch(
      #
      #  Multiplying a matrix by a diagonal matrix can be replicated by multiplying each
      #  column by the diagonal's corresponding entry and thus avoids creating and storing
      #  the diagonal matrix resulting in much more efficient memory usage. Grouping the 
      #  last few terms results in fewer calculations and minor speed savings. 
      #
      solve(t(X * V) %*% X) %*% (t(X) %*% (V * y)), 
          
      error=function(e) { 
        #  
        #  Catch errors relating to non-conformable dimensions by examining the error 
        #  message and then throw our own error with a more helpful message.
        #
        if (e$message=='non-conformable arguments'| grepl('dims', e$message)) {
          stop(paste0('Dimensions of initial arguments are non-conformable:',
                      '\n  -  y must be an nx1 vector, currently ', NROW(y), 'x', NCOL(y),  
                      '\n  -  X must be an nxp matrix, currently ', NROW(X), 'x', NCOL(X),
                      '\n  -  V must be an nx1 vector, currently ', NROW(V), 'x', NCOL(V)))
          
        #  If an unrelated error arises allow R to proceed as it would otherwise. 
        } else stop(e)}),    
    
      #  Use invokeRestart function to suppress warnings in the case we throw an error.
      warning=function(w) {invokeRestart('muffleWarning')})
  #
  #  Calculate siqsq, the weighted estimator for the standard error. We need to force
  #  sigsq to be 'numeric' in order to compute the coefficient estimators' covariance
  #  matrix later on.
  #
  sigsq <- (t(y - X %*% betahat) * V) %*% (y - X %*% betahat) / (n-p)
  sigsq <- as.numeric(sigsq)
  #
  #  Calculate betahat.cov, the coefficient estimators' covariance matrix, and then 
  #  define sebeta, the standard error for each estimator, as its diagonal.
  #
  betahat.cov <- solve(t(X * V) %*% X) * sigsq
  sebeta <- as.matrix(diag(betahat.cov))
  #
  #  Finally, return a named list of the results.  
  #
  list(betahat=betahat, sebeta=sebeta, sigsq=sigsq)
}