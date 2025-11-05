### Methods for icfit objects

setMethod("show", signature(object='icfit'), function(object) {
  message("icfit object\n\n")
  message("Model with best ",object$criterion,": ",paste(object$topmodel, collapse=' '),"\n\n")
  message("Use summary(), coef() and predict() to get inference for the top model\n")
  message("Use coef(object$msfit) and predict(object$msfit) to get BMA estimates and predictions\n")
}
)

confint.icfit <- function(object, ...) {
    confint(object$topmodel.fit)
}

coef.icfit <- function(object,...) {
    coef(object$topmodel.fit)
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
    new("icfit",ans)
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

findmodels_fast <- function(y, x, data, smoothterms, nknots=9, groups, constraints, enumerate, includevars, maxvars, niter=5000, family='normal', priorCoef, priorGroup, priorModel=modelbbprior(1,1), priorConstraints, priorVar=igprior(.01,.01), priorSkew=momprior(tau=0.348), neighbours, phi, deltaini, adj.overdisp='intercept', hess='asymp', optimMethod, optim_maxit, initpar='none', B=10^5, XtXprecomp, verbose=TRUE, fastmethod) {
  if (!fastmethod %in% c("L0Learn","L1","CDA","all")) stop("fastmethod must be 'L0Learn', 'L1', 'CDA' or 'all'")
  #Determine what model search methods to use
  L0Learn_available <- ifelse(family %in% c("normal","binomial"), TRUE, FALSE)
  L1_available <- ifelse(family %in% c("normal","twopiecenormal","laplace","twopiecelaplace","auto","binomial","binomial logit","poisson","poisson log","negbinomial","multinomial","cox","mnormal"), TRUE, FALSE)
  CDA_available <- ifelse(family %in% c("normal","twopiecenormal","laplace","twopiecelaplace","auto","binomial","binomial logit","poisson","poisson log"), TRUE, FALSE)
  use_L0Learn <- ifelse((fastmethod %in% c("L0Learn","all")) & L0Learn_available, TRUE, FALSE)
  use_L1 <- ifelse((fastmethod %in% c("L1","all")) & L1_available, TRUE, FALSE)
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
  models <- rbind(models_L0Learn, models_L1, models_CDA)  
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


findmodels_lasso <- function(y, x, data, smoothterms, nknots=9, groups, constraints, enumerate, includevars, maxvars, niter=5000, family='normal', priorCoef, priorGroup, priorModel=modelbbprior(1,1), priorConstraints, priorVar=igprior(.01,.01), priorSkew=momprior(tau=0.348), neighbours, phi, deltaini, adj.overdisp='intercept', hess='asymp', optimMethod, optim_maxit, initpar='none', B=10^5, XtXprecomp, verbose=TRUE) {
  # Build design matrix
  tmp <- modelSelection:::formatInputdata(y=y,x=x,data=data,smoothterms=smoothterms,nknots=nknots,family=family)
  x <- tmp$x; y <- tmp$y
  # Call glmnet
  family_glmnet <- family2glmnet(family)
  fit <- glmnet::glmnet(x=x, y=y, family=family_glmnet, alpha=1, intercept=FALSE, standardize=TRUE)
  models <- unique(t(as.matrix(coef(fit) != 0)))
  return(models)
}



bestBIC_fast <- function(..., fastmethod) {
  args <- list(...)
  args$models <- findmodels_fast(..., priorCoef=bic(), priorModel=modelunifprior(), fastmethod=fastmethod)
  ans <- do.call(bestBIC, args)
  #ans <- bestEBIC(..., models=models)
  return(ans)  
}

bestAIC_fast <- function(..., fastmethod) {
  args <- list(...)
  args$models <- findmodels_fast(..., priorCoef=aic(), priorModel=modelunifprior(), fastmethod=fastmethod)
  ans <- do.call(bestAIC, args)
  return(ans)  
}

bestEBIC_fast <- function(..., fastmethod) {
  args <- list(...)
  args$models <- findmodels_fast(..., priorCoef=bic(), priorModel=modelbbprior(), fastmethod=fastmethod)
  ans <- do.call(bestEBIC, args)
  return(ans)  
}

bestIC_fast <- function(..., penalty, fastmethod) {
  args <- list(...)
  args$models <- findmodels_fast(..., priorCoef=ic(penalty), priorModel=modelunifprior(), fastmethod=fastmethod)
  ans <- do.call(bestIC, args, penalty)
  return(ans)  
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



