### Methods for icfit objects

setMethod("show", signature(object='icfit'), function(object) {
  cat("icfit object\n\n")
  cat("Model with best ",object$criterion,": ",paste(object$topmodel, collapse=' '),"\n\n")
  cat("Use summary(), coef() and predict() to get inference for the top model\n")
  cat("Use coef(object$msfit) and predict(object$msfit) to get BMA estimates and predictions\n")
}
)

confint.icfit <- function(object, ...) {
    confint(object$topmodel.fit)
}

#Return non-zero coefficient estimates
coef.icfit <- function(object,...) {
    coef(object$topmodel.fit)
}

#Return coefficient estimates for all covariates, including the non-zeroes
coefall <- function(object) {  UseMethod("coefall") }

coefall.icfit <- function(object) {
  b <- rep(0, length(object$varnames))
  names(b) <- object$varnames
  b[object$topmodel] <- coef(object$topmodel.fit)    
  return(b)
}


predict.icfit <- function(object, ...) {
    predict(object$topmodel.fit, ...)
}

summary.icfit <- function(object, ...) {
    summary(object$topmodel.fit, ...)
}




## FIND THE MODEL ATTAINING THE BEST VALUE OF AN INFORMATION CRITERION

checkargs_IC <- function(...) {
    params= eval(quote(list(...)))
    forbiddenpars= c('priorCoef','priorModel','center','scale')
    if (any(forbiddenpars %in% names(params))) stop(paste("Arguments",paste(forbiddenpars,collapse=", "),"have set values so they cannot be passed on to modelSelection"))
}


topmodelnames <- function(ans) {
    topvarids= as.numeric(strsplit(ans$models$modelid[1], ',')[[1]])
    ans= list(topvarids= topvarids, varnames= ans$varnames[topvarids])
    return(ans)
}

family2glm <- function(family) {
    if (family=="normal") {
        ans= gaussian()
    } else if (family=="binomial") {
        ans= binomial()
    } else if (family=="poisson") {
        ans= poisson()
    } else stop("Only the Normal, Binomial and Poisson families are currently implemented")
    return(ans)
}


#Extract info criteria from an msfit and return icfit object
extractmsIC <- function(ms, getICfun) {
    ans= vector("list",5)
    names(ans)= c('topmodel','topmodel.fit','models','varnames','msfit')
    ans$models= getICfun(ms)
    if (!is.null(colnames(ms$xstd))) {
        ans$varnames= colnames(ms$xstd)
    } else {
        ans$varnames= paste('x[,',1:ncol(ms$xstd),']',sep='')
    }
    tm= topmodelnames(ans)
    ans$topmodel= tm$varnames
    ans$msfit= ms
    data= data.frame(y=ms$ystd, ms$xstd[,tm$topvarids])
    names(data)[-1]= ans$varnames[tm$topvarids]
    f= formula(y ~ -1 + .)
    ans$topmodel.fit= glm(f , data=data, family=family2glm(ms$family))
    new("icfit", ans)
}

#Combine two icfit objects into a single one
combineICfit <- function(fit1, fit2) {
    if (any(fit1$varnames != fit2$varnames)) stop("Cannot combine two icfit objects with different varnames")
    ans= vector("list",5)
    names(ans)= c('topmodel','topmodel.fit','models','varnames','msfit')
    ans$varnames = fit1$varnames
    ans$models <- rbind(fit1$models, fit2$models)
    ans$models <- dplyr::arrange(ans$models, ans$models[[2]])  #sort increasingly by value of the information criterion
    if (fit1$models[1,2] < fit2$models[1,2]) {
        ans$topmodel = fit1$topmodel
        ans$topmodel.fit = fit1$topmodel.fit
    } else {
        ans$topmodel = fit2$topmodel
        ans$topmodel.fit = fit2$topmodel.fit
    }
    new("icfit", ans)
}



bestBIC <- function(...) {
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=bic(), priorModel=modelunifprior(), center=FALSE, scale=FALSE)
    ans= extractmsIC(ms, getBIC)
    ans$criterion= 'BIC'
    return(ans)
}

bestAIC <- function(...) {
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=aic(), priorModel=modelunifprior(), center=FALSE, scale=FALSE)
    ans= extractmsIC(ms, getAIC)
    ans$criterion= 'AIC'
    return(ans)
}

bestEBIC <- function(...) {
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=bic(), priorModel=modelbbprior(), center=FALSE, scale=FALSE)
    ans= extractmsIC(ms, getEBIC)
    ans$criterion= 'EBIC'
    return(ans)
}

bestIC <- function(..., penalty) {
    if (missing(penalty)) stop("penalty must be specified. Alternatively consider using bestBIC(), bestEBIC() or bestAIC()")
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=ic(penalty), priorModel=modelunifprior(), center=FALSE, scale=FALSE)
    ans= extractmsIC(ms, getIC)
    ans$criterion= paste('GIC (penalty=',penalty,')',collapse='')
    return(ans)
}

bestBIC_fast <- function(..., fastmethod='all') {
  args <- list(...)
  args$models <- findmodels_fast(..., priorCoef=bic(), priorModel=modelunifprior(), fastmethod=fastmethod)
  ans <- do.call(bestBIC, args)
  if (fastmethod %in% c('adaptiveL1','all')) {
    newmodels <- findmodels_alasso(ans, ...)
    if (!is.null(newmodels)) {
      args$models <- newmodels
      ansnew <- do.call(bestBIC, args)
      ans <- combineICfit(ans, ansnew)
    }
  } 
  return(ans)  
}

bestAIC_fast <- function(..., fastmethod='all') {
  args <- list(...)
  args$models <- findmodels_fast(..., priorCoef=aic(), priorModel=modelunifprior(), fastmethod=fastmethod)
  ans <- do.call(bestAIC, args)
  if (fastmethod %in% c('adaptiveL1','all')) {
    newmodels <- findmodels_alasso(ans, ...)
    if (!is.null(newmodels)) {
      args$models <- newmodels
      ansnew <- do.call(bestAIC, args)
      ans <- combineICfit(ans, ansnew)
    }
  } 
  return(ans)  
}

bestEBIC_fast <- function(..., fastmethod='all') {
  args <- list(...)
  args$models <- findmodels_fast(..., priorCoef=bic(), priorModel=modelbbprior(), fastmethod=fastmethod)
  ans <- do.call(bestEBIC, args)
  if (fastmethod %in% c('adaptiveL1','all')) {
    newmodels <- findmodels_alasso(ans, ...)
    if (!is.null(newmodels)) {
      args$models <- newmodels
      ansnew <- do.call(bestEBIC, args)
      ans <- combineICfit(ans, ansnew)
    }
  } 
  return(ans)  
}

bestIC_fast <- function(..., penalty, fastmethod='all') {
  args <- list(...)
  args$models <- findmodels_fast(..., priorCoef=ic(penalty), priorModel=modelunifprior(), fastmethod=fastmethod)
  ans <- do.call(bestIC, args, penalty)
  if (fastmethod %in% c('adaptiveL1','all')) {
    newmodels <- findmodels_alasso(ans, ...)
    if (!is.null(newmodels)) {
      args$models <- newmodels
      ansnew <- do.call(bestIC, args, penalty)
      ans <- combineICfit(ans, ansnew)
    }
  } 
  return(ans)  
}


findmodels_fast <- function(y, x, data, smoothterms, nknots=9, groups, constraints, enumerate, includevars, maxvars, niter=5000, family='normal', priorCoef, priorGroup, priorModel=modelbbprior(1,1), priorConstraints, priorVar=igprior(.01,.01), priorSkew=momprior(tau=0.348), neighbours, phi, deltaini, adj.overdisp='intercept', hess='asymp', optimMethod, optim_maxit, initpar='none', B=10^5, XtXprecomp, verbose=TRUE, fastmethod) {
  if (!fastmethod %in% c("L0Learn","L1","adaptiveL1","CDA","all")) stop("fastmethod must be 'L0Learn', 'L1', 'adaptiveL1', 'CDA' or 'all'")
  #Determine what model search methods to use
  L0Learn_available <- ifelse(family %in% c("normal","binomial"), TRUE, FALSE)
  L1_available <- ifelse(family %in% c("normal","twopiecenormal","laplace","twopiecelaplace","auto","binomial","binomial logit","poisson","poisson log","negbinomial","multinomial","cox","mnormal"), TRUE, FALSE)
  CDA_available <- ifelse(family %in% c("normal","twopiecenormal","laplace","twopiecelaplace","auto","binomial","binomial logit","poisson","poisson log"), TRUE, FALSE)
  use_L0Learn <- ifelse((fastmethod %in% c("L0Learn","all")) & L0Learn_available, TRUE, FALSE)
  use_L1 <- ifelse((fastmethod %in% c("L1","adaptiveL1","all")) & L1_available, TRUE, FALSE)
  use_CDA <- ifelse((fastmethod %in% c("CDA","all")) & CDA_available, TRUE, FALSE)
  if (!use_L0Learn & !use_L1 & !use_CDA) stop("The specified fastmethod is unavailable for the specified family")
  #Run L0Learn
  if (use_L0Learn) {
    loss <- ifelse(family == "normal", "SquaredError", "Logistic")
    # Build design matrix
    tmp <- formatInputdata(y=y,x=x,data=data,smoothterms=smoothterms,nknots=nknots,family=family)
    x <- tmp$x; y <- tmp$y
    # Call L0Learn
    if (missing(includevars)) includevars <- rep(FALSE, ncol(x))
    if (missing(maxvars)) maxvars <- set_ms_maxvars(n=nrow(x), p=ncol(x), priorCoef=bicprior(), family=family, includevars=includevars)
    fit <- L0Learn.fit(x=x, y=y, penalty="L0", maxSuppSize=maxvars, loss=loss, algorithm="CDPSI")
    models <- unique(t(as.matrix(coef(fit) != 0)))
    if (ncol(models) > ncol(x)) models <- models[,-1] #remove intercept
    models_L0Learn <- models
  } else {
    models_L0Learn <- NULL
  }
  #Run L1
  if (use_L1) {
    models_L1 <- findmodels_lasso(y=y, x=x, data=data, smoothterms, nknots=nknots, groups=groups, constraints=constraints, enumerate=enumerate, includevars=includevars, maxvars=maxvars, niter=niter, family=family, priorCoef=priorCoef, priorGroup=priorGroup, priorModel=priorModel, priorConstraints=priorConstraints, priorVar=priorVar, priorSkew=priorSkew, neighbours=neighbours, phi=phi, deltaini=deltaini, adj.overdisp=adj.overdisp, hess=hess, optimMethod=optimMethod, optim_maxit=optim_maxit, initpar=initpar, B=B, XtXprecomp=XtXprecomp, verbose=verbose)
  } else {
    models_L1 <- NULL
  }
  #Run CDA
  if (use_CDA) {
    fit <- modelSelection(y=y, x=x, data=data, initSearch='CDA', burnin=0, niter=1, center=FALSE, scale=FALSE, smoothterms=smoothterms, nknots=nknots, groups=groups, constraints=constraints, enumerate=enumerate, includevars=includevars, maxvars=maxvars, priorCoef=priorCoef, priorGroup=priorGroup, priorModel=priorModel, priorConstraints=priorConstraints, priorVar=priorVar, priorSkew=priorSkew, neighbours=neighbours, phi=phi, deltaini=deltaini, adj.overdisp=adj.overdisp, hess=hess, optimMethod=optimMethod, optim_maxit=optim_maxit, initpar=initpar, XtXprecomp=XtXprecomp, verbose=verbose)
    models_CDA <- matrix(fit$postMode == 1, nrow=1)
  } else {
    models_CDA <- NULL
  }
  models <- unique(rbind(models_L0Learn, models_L1, models_CDA))
  return(models)
}


# Map GLM family from the format required by modelSelection to that used by glmnet
# - Families twopiecenormal, laplace and twopiecelaplace are mapped to normal, since these are unavailable in glmnet
# - Family negbinomial is mapped to Poisson, since the negative binomial is unavailable in glmnet
family2glmnet <- function(family) {
    if (family %in% c("normal","twopiecenormal","laplace","twopiecelaplace","auto")) {
        family_glmnet = "gaussian"
    } else if (family %in% c("binomial","binomial logit")) {
        family_glmnet = "binomial"
    } else if (family %in% c("poisson", "poisson log", "negbinomial")) {
        family_glmnet = "poisson"
    } else if (family == "multinomial") {
        family_glmnet = "multinomial"
    } else if (family == "cox") {
        family_glmnet = "cox"
    } else if (family == "mgaussian") {
        family_glmnet = "mgaussian"
    } else {
        stop("This family is unavailable in glmnet")
    }
    return(family_glmnet)
}


# Return all models visited by the LASSO regularization path. If bini is not missing, adaptive LASSO is used (only on variables such that bini != 0)
findmodels_lasso <- function(bini, y, x, data, smoothterms, nknots=9, groups, constraints, enumerate, includevars, maxvars, niter=5000, family='normal', priorCoef, priorGroup, priorModel=modelbbprior(1,1), priorConstraints, priorVar=igprior(.01,.01), priorSkew=momprior(tau=0.348), neighbours, phi, deltaini, adj.overdisp='intercept', hess='asymp', optimMethod, optim_maxit, initpar='none', B=10^5, XtXprecomp, verbose=TRUE) {
  # Build design matrix
  tmp <- formatInputdata(y=y,x=x,data=data,smoothterms=smoothterms,nknots=nknots,family=family)
  x <- tmp$x; y <- tmp$y
  family_glmnet <- family2glmnet(family)
  if (missing(bini)) { # LASSO
    fit <- glmnet::glmnet(x=x, y=y, family=family_glmnet, alpha=1, intercept=FALSE, standardize=TRUE)
    models <- unique(t(as.matrix(coef(fit)[-1,,drop=FALSE] != 0)))
  } else {             # Adaptive LASSO based on initial parameter estimate bini
    sel <- (bini != 0)
    xstd <- scale(x[,sel,drop=FALSE])
    xstd <- t(t(xstd) * abs(bini[sel]))  
    fit <- glmnet::glmnet(x=xstd, y=y, family=family_glmnet, alpha=1, intercept=FALSE, standardize=FALSE)
    msel <- unique(t(as.matrix(coef(fit)[-1,,drop=FALSE] != 0)))
    models <- matrix(FALSE, nrow=nrow(msel), ncol=ncol(x))
    models[,sel] <- msel
  }
  return(models)
}


# Given icfit object, search for new models using adaptive LASSO and add them to the icfit object
findmodels_alasso <- function(ans, ...) {
  newmodels <- NULL
  if (length(coef(ans)) > 1) {
    bini <- coefall(ans)
    newmodels <- findmodels_lasso(bini, ...)
    newmodelid <- apply(newmodels, 1, function(z) paste(which(z), collapse = ","))
    # Keep only new models that aren't already in args$models
    isnew <- !(newmodelid %in% ans$models$modelid)
    if (any(isnew)) {
      newmodels <- newmodels[isnew,,drop=FALSE]
    }
  }
  return(newmodels)
}


## EXTRACT INFORMATION CRITERIA FROM AN msfit object

setMethod("getAIC", signature(object='msfit'), function(object) {
    ans= getBIC(object)
    names(ans)[names(ans)=='bic']= 'aic'
    return(ans)
}
)


setMethod("getBIC", signature(object='msfit'), function(object) {
    pc= object$priors$priorCoef
    pm= object$priors$priorModel
    if ((pc@priorDistr != 'bic') || (pm@priorDistr != 'uniform')) stop("To obtain BIC you should set priorCoef=bic() and priorModel=modelunifprior() when calling modelSelection")
    ans= getIC(object)
    names(ans)[names(ans)=='ic']= 'bic'
    return(ans)
}
)

setMethod("getIC", signature(object='msfit'), function(object) {
    ic = -2 * ( object$postProb + object$p * log(2) )
    ans= dplyr::tibble(modelid=object$modelid, ic=ic)
    ans= ans[order(ans$ic),]
    return(ans)
}
)


setMethod("getEBIC", signature(object='msfit'), function(object) {
    pc= object$priors$priorCoef
    pm= object$priors$priorModel
    isbbprior= (pm@priorDistr == 'binomial') && all(pm@priorPars == c(1,1))
    if ((pc@priorDistr != 'bic') || !isbbprior) stop("To obtain BIC you should set priorCoef=bic() and priorModel=modelbbprior() when calling modelSelection")
    ebic = -2  * ( object$postProb + log(object$p+1) )
    ans= dplyr::tibble(modelid=object$modelid, ebic=ebic)
    ans= ans[order(ans$ebic),]
    return(ans)
}
)



