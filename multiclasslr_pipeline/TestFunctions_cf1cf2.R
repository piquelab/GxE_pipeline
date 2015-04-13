## given a design matrix & class counts, returns a vector of zero weights
InitialWeights <- function(design.mat, n.classes)
{
  n.params.ni <- dim(design.mat)[2]
  initweights.ni <- rep(0, n.params.ni*n.classes)
  initweights.ni
}

## given a design matrix & class labels, returns a vector of complete variable names
VarNames <- function(design.mat, configurations)
{
  indvars <- colnames(design.mat)
  if(is.null(indvars)){stop("design matrix must have column names")}
  var.names <- matrix(sapply(configurations, function(config){t(paste0(indvars, config))}), ncol=1)[, 1]
  var.names
}
 
## multinomial data will be transformed for ggplot2
CreatePlotDat <- function(multiDat, sig.parms=FALSE)
{
  if(!("results" %in% names(multiDat))){stop("\"results\" absent from input data")}
  if(!("dat.all" %in% names(multiDat$results))){stop("\"dat.all\" absent from results data")}

  if(sig.parms==TRUE){
    dat <- multiDat$results$dat.sig
  } else {
    dat <- multiDat$results$dat.all
  }
  
  var.names <-  rownames(dat)
  all.data <- data.frame(var.names, dat, stringsAsFactors=FALSE)
  intercept_exclude <- !(all.data$var.names %in% c('intercept_cf1', 'intercept_cf2', 'intercept_cf3'))
  feature <- sapply(all.data$var.names, function(id){strsplit(id, '_')[[1]][1]}, USE.NAMES=FALSE)
  config <- as.factor(sapply(all.data$var.names, function(id){strsplit(id, '_')[[1]][2]}, USE.NAMES=FALSE))
  plot_dat <- data.frame(feature, config, parms=all.data$parms, se=all.data$parms.se, stringsAsFactors=FALSE)[intercept_exclude, ] 
  plot_dat
}

## clojure that bundles a baselined design matrix into a function that returns a logistic or linear regression
## depending on the input
PrepLinearModel <- function(X.baseline)
{
  MyLinearModel <- function(Y, X=X.baseline)
  {
    ##Y <- p.1_lm
    ##X <- X.baseline
    exclude <- is.na(Y)
    Y <- Y[!exclude]
    X <- X[!exclude, ] 
    dat <- data.frame(Y, X)
    if(sum(dat$Y < 0) == 0){ 
      mod_logit <- glm(Y ~ ., data=dat, family=binomial(logit))
      mod_logit
    } else {
      mod_lm <- lm(Y ~ ., data=dat)
      mod_lm
    }
  }
  MyLinearModel
}

## log-likelihood closure, bundling the intercept & number of classes into a likelihood function
## for optim
MultiNegLlk.ni <- function(Y, X, intercepts, n.classes)
{
  Y <- as.matrix(Y)
  X <- as.matrix(X)  
  MyMultiNegLlk <- function(w)
  {
    # add an intercept term & intercept weights  
    X.int <- cbind(rep(1, dim(X)[1]), X)
    #w.int <- c(intercepts[1], w[1:37], intercepts[2], w[38:74], intercepts[3], w[75:111])
    parm.classes <- n.classes-1
    n.all.parms <- length(w)
    if((n.all.parms %% parm.classes) != 0){stop("parms vector length is not a multiple of n.classes-1")}
    n.parms <- n.all.parms / parm.classes
    w.int <- c(intercepts[1], w[1:n.parms], intercepts[2], w[(n.parms+1):(2*n.parms)], intercepts[3], w[(2*n.parms+1):(3*n.parms)])

    Y_hat <- YHatMulti(w.int, X.int)  
    llk <- sum(diag(t(Y) %*% log(Y_hat)))  # trace is invariant to cyclic permutations
    -llk
  }  
  MyMultiNegLlk
}

## constrain parameters for cf1 & cf2 log-likelihood closure, bundling the intercept & number
## of classes into a likelihood function for optim
MultiNegLlk.ni.cf1cf2 <- function(Y, X, intercepts, n.classes)
{
  Y <- as.matrix(Y)
  X <- as.matrix(X)  
  MyMultiNegLlk <- function(w)
  {
      fun.start <- proc.time()
    # add an intercept term & intercept weights  
    X.int <- cbind(rep(1, dim(X)[1]), X)
    #w.int <- c(intercepts[1], w[1:37], intercepts[2], w[38:74], intercepts[3], w[75:111])
    #parm.classes <- n.classes-1
    ## remove the baseline class & combine cf1 & cf2
    parm.classes <- n.classes-2
    n.all.parms <- length(w)
    #if((n.all.parms %% parm.classes) != 0){stop("parms vector length is not a multiple of n.classes-1")}
    if((n.all.parms %% parm.classes) != 0){stop("parms vector length is not a multiple of n.classes-2")}
    n.parms <- n.all.parms / parm.classes
    w.cf1cf2 <- c(intercepts[1], w[1:n.parms])
    w.int <- c(rep(w.cf1cf2, 2), intercepts[2], w[(n.parms+1):(2*n.parms)])
      yhat.start <- proc.time()
    Y_hat <- YHatMulti(w.int, X.int)  
      yhat.end <- proc.time()
      yhat.cumtime <<- yhat.cumtime + (yhat.end - yhat.start)
    ## here just do the operations on the diagnoal
    llk <- sum(diag(t(Y) %*% log(Y_hat)))  # trace is invariant to cyclic permutations
      fun.end <- proc.time()
      fun.cumtime <<- fun.cumtime + (fun.end - fun.start)
    #cat('cumtime: ', cumtime, '\n')
    -llk
  }
  MyMultiNegLlk
}

## optimization for a multinomial model. intercept fit from an independent
## run, avoiding the use of a baseline class
MyMultinomial.ni <- function(intercepts, parms, X, Y, var.names, hessian=TRUE, maxit=1E5, cf1cf2=FALSE)
{
  ## check argument validity ##
  ## TODO ##
  n.classes <- dim(Y)[2]
  cat('=== parameters are valid ===', '\n\n')

  ## create gradient closure with intercepts
  if(cf1cf2==FALSE){
    MyMultiNegLlk <- MultiNegLlk.ni(Y, X, intercepts, n.classes)
  } else {
    MyMultiNegLlk <- MultiNegLlk.ni.cf1cf2(Y, X, intercepts, n.classes)
  }
    
  ## model fitting
  cat('=== begin optimization ===', '\n\n')
    run.start <- proc.time()
    fun.cumtime <- run.start * 0
    yhat.cumtime <- fun.cumtime
  model <- optim(parms, fn=MyMultiNegLlk, method="BFGS", hessian=TRUE, control=list(maxit=maxit));
    run.end <- proc.time()
  cat('Run time: ', run.end-run.start, '\n\n');
  cat('Function time: ', fun.cumtime, '\n');
  cat('Yhat time: ', yhat.cumtime, '\n');
  ## results: model fit and parameter estimates
  cat('=== processing ouput ===', '\n\n')
  nllk <- model$value 
  converged <- model$convergence
  hessian <- model$hessian
  model_res <- MungeResults(model, var.names, classes=(n.classes-1))

  retval <- list(nllk=nllk, converged=converged, hessian=hessian, results=model_res)
  retval
}

#########################################################################################
## adapt the multinomial function to handle intercepts on the constrained model,
## combine cf1 & cf2
#########################################################################################
MyMultinomial <- function (parms, X, Y, var.names, hessian = TRUE, maxit = 1e+05, io.cf1cf2=FALSE) 
{
    cat("=== parameters are valid ===", "\n\n")
    cat("=== begin optimization ===", "\n\n")

    if(io.cf1cf2==FALSE){
      MyLlk <- MultiNegLlk
    } else {
      MyLlk <- MultiNegLlk.io.cf1cf2
    }
    
    st <- proc.time()
    model <- optim(parms,
                   fn = MyLlk,
                   X = as.matrix(X), 
                   Y = as.matrix(Y),
                   method = "BFGS",
                   hessian = TRUE,
                   control = list(maxit = maxit))
    
    en <- proc.time()

    cat("Run time: ", en - st, "\n\n")
    cat("=== processing ouput ===", "\n\n")
    nllk <- model$value
    converged <- model$convergence
    hessian <- model$hessian
    model_res <- MungeResults(model, var.names, classes = 4)
    retval <- list(nllk = nllk, converged = converged, hessian = hessian, 
        results = model_res)
    retval
}

MultiNegLlk.io.cf1cf2 <- function (w, X, Y) 
{  
    w.io.cf1cf2 <- c(rep(w[1], 2), w[2])  
    Y_hat <- YHatMulti(w.io.cf1cf2, X)
    llk <- sum(diag(t(Y) %*% log(Y_hat)))
    -llk
}

##         ##
## THE END ##
##         ##
