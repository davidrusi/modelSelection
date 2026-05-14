###
### ggm.R
###

### Methods for msfit_ggm objects

plot.msfit_ggm= function(x, y, ...) {
  postSample= x$postSample[,x$indexes[1,] != x$indexes[2,],drop=FALSE]
  margppcum= apply(postSample !=0, 2, cumsum) / seq_len(nrow(postSample))
  plot(margppcum[,1], type='l', ylim=c(0,1), xlab='Iteration', ylab='Marginal posterior inclusion probabilities')
  if (ncol(margppcum)>1) for (i in 2:min(5000,ncol(margppcum))) lines(margppcum[,i])
}


setMethod("show", signature(object='msfit_ggm'), function(object) {
  message('Gaussian graphical model (msfit_ggm object) with ',object$p,' variables\n')
  message("Use coef() to get BMA estimates, posterior intervals and posterior marginal prob of entries being non-zero\n")
  message("use postProb() for posterior model probabilities")
}
)

coef.msfit_ggm <- function(object,...) {
  m= as.vector(object$postmean) #use Rao-Blackwellized estimator
  #m= Matrix::colMeans(object$postSample) #use MCMC estimator
  ci= sparseMatrixStats::colQuantiles(object$postSample, prob=c(0.025,0.975))
  ans= cbind(t(object$indexes), m, ci, object$margpp) #use Rao-Blackwellized edge inclusion probabilities
  colnames(ans)[-1:-2]= c('estimate','2.5%','97.5%','margpp')
  return(ans)
}

icov <- function(fit, threshold) {
  if (!inherits(fit, 'msfit_ggm')) stop("Argument fit must be of class msfit_ggm")
  m= as.vector(fit$postmean) #use Rao-Blackwellized estimator
  #m= Matrix::colMeans(fit$postSample) #use MCMC estimator
  if (!missing(threshold)) {
      sel= (fit$margpp >= threshold) #use Rao-Blackwellized edge inclusion probabilities
      #sel= (Matrix::colMeans(fit$postSample != 0) >= threshold) #use MCMC edge inclusion probabilities
      ans= Matrix::sparseMatrix(i= fit$indexes[1,sel], j= fit$indexes[2,sel], x=m[sel], dims=c(fit$p,fit$p), symmetric=TRUE)
  } else {
      ans= matrix(nrow=fit$p,ncol=fit$p)
      ans[upper.tri(ans,diag=TRUE)]= m
      ans[lower.tri(ans)]= t(ans)[lower.tri(ans)]
  }
  return(ans)
}


proportion_visited_models= function(postSample) {
  modelpp <- apply(postSample != 0, 1, function(z) paste(which(z),collapse=','))
  modelpp <- table(modelpp)/length(modelpp)
  modelpp <- data.frame(modelid=names(modelpp), pp=as.numeric(modelpp))
  return(modelpp)
}

setMethod("postProb", signature(object='msfit_ggm'), function(object, nmax, method='norm') {
if (!is.null(object$models)) {
  ans= object$models
} else {
  param_ids= paste('(',object$indexes[1,], ',', object$indexes[2,], ')',sep='')
  npar= ncol(object$indexes)
  #Compute posterior probabilities
  modelpp = proportion_visited_models(object$postSample)
  #Compute model identifiers
  tmp= strsplit(modelpp$modelid, ',')
  modelid= t(sapply(tmp, function(z) { ans= rep(FALSE,npar); ans[as.integer(z)]= TRUE; return(ans) }))
  colnames(modelid)= param_ids
  modelid= Matrix::Matrix(modelid, sparse=TRUE)
  #Sort models in decreasing probability
  o= order(modelpp$pp,decreasing=TRUE)
  modelpp= modelpp[o,]
  modelid= modelid[o,]
  #Return nmax models    
  if (!missing(nmax)) {
      modelpp <- modelpp[1:nmax,]
      modelid = modelid[1:nmax,]
  }
  #Return output
  ans= list(modelid= modelid, pp=modelpp$pp)
}
return(ans)
}
)



### Model selection routines

modelSelectionGGM= function(y, priorCoef=normalidprior(tau=1), priorModel=modelbinomprior(1/ncol(y)), priorDiag=exponentialprior(lambda=1), center=TRUE, scale=TRUE, global_proposal= "none", prob_global=0.5, tempering=0.5, truncratio= 100, save_proposal=FALSE, niter=10^3, burnin= round(niter/10), updates_per_iter= ceiling(sqrt(ncol(y))), updates_per_column= 10, sampler='Gibbs', pbirth=0.75, pdeath=0.5*(1-pbirth), bounds_LIT, Omegaini='glasso-ebic', verbose=TRUE) {
  #Check input args
  if (!is.matrix(y)) y = as.matrix(y)
  p= ncol(y);
  if (p <=1) stop("y must have at least 2 columns")
  if (!is.numeric(y)) stop("y must be numeric")
  if (!(sampler %in% c('Gibbs','birthdeath','LIT'))) stop("sampler must be 'Gibbs', 'birthdeath', or 'LIT'")
  if (tempering < 0) stop("tempering cannot be negative")
  y = scale(y, center=center, scale=scale)
    
  #Format prior parameters
  prCoef= formatmsPriorsMarg(priorCoef=priorCoef, priorVar=priorDiag)
  prCoef= as.list(c(priorlabel=prCoef$priorCoef@priorDistr, prCoef[c('prior','tau','lambda')]))
  prModel= as.list(c(priorlabel=priorModel@priorDistr, priorPars= priorModel@priorPars))
  prModel[grep("priorPars", names(prModel))] <- as.double(prModel[grep("priorPars", names(prModel))])
    
  #Format posterior sampler parameters
  if (missing(bounds_LIT)) bounds_LIT= log(c("lbound_death"=1/p, "ubound_death"= 1, "lbound_birth"=1/p, "ubound_birth"=p))
  samplerPars= format_GGM_samplerPars(sampler, p=p, niter=niter, burnin=burnin, updates_per_iter=updates_per_iter, updates_per_column = updates_per_column, pbirth=pbirth, pdeath=pdeath, prob_global=prob_global, tempering=tempering, truncratio=truncratio, global_proposal=global_proposal, bounds_LIT=bounds_LIT, verbose=verbose)
    
  #Initial value for sampler
  if (verbose) message(" Obtaining initial parameter estimate...")
  Omegaini= initialEstimateGGM(y, Omegaini)
  if (verbose) message(" Done\n")
    
  #Call C++ function
  proposal= proposaldensity= NULL
  if (global_proposal == 'none') {
    ans= modelSelectionGGMC(y, prCoef, prModel, samplerPars, Omegaini)
    postSample= Matrix::t(ans$postSample)
    postmean = ans$postmean
    margpp= ans$margpp
    prop_accept= ans$prop_accept
  } else {
    ans= modelSelectionGGM_globalC(y, prCoef, prModel, samplerPars, Omegaini)
    postSample= Matrix::t(ans[[1]])
    postmean= ans[[2]]
    margpp= ans[[3]]
    prop_accept= ans[[4]]
    if (save_proposal) {
      proposal= lapply(ans[5:(5+p-1)], Matrix::t)
      for (i in seq_len(length(proposal))) proposal[[i]]@x = as.double(proposal[[i]]@x)
      proposaldensity= Matrix::t(ans[[length(ans)]])
    }
  }
    

  #Return output
  priors= list(priorCoef=priorCoef, priorModel=priorModel, priorDiag=priorDiag)

  A= diag(p)
  indexes= rbind(row(A)[upper.tri(row(A),diag=TRUE)], col(A)[upper.tri(row(A),diag=TRUE)])
  rownames(indexes)= c('row','column')

  ans= list(postSample=postSample, prop_accept=prop_accept, proposal=proposal, proposaldensity=proposaldensity, postmean=postmean, margpp=margpp, priors=priors, p=p, indexes=indexes, samplerPars=samplerPars, global_proposal=global_proposal)

  new("msfit_ggm",ans)
}



#Format posterior sampling parameters to pass onto C++
format_GGM_samplerPars= function(sampler, p, niter, burnin, updates_per_iter, updates_per_column, pbirth, pdeath, prob_global, tempering, truncratio, global_proposal, bounds_LIT, verbose) {
  if (missing(updates_per_column)) {
      updates_per_column= as.integer(min(p, max(log(p), 10)))
  } else {
      updates_per_column= as.integer(updates_per_column)
  }
  if (!(global_proposal %in% c('none','regression','in-sample'))) stop("global_proposal must be 'none', 'regression' or 'in-sample'")
  if ((updates_per_iter < 1) || (updates_per_iter > p)) stop("updates_per_iter must be between 1 and ncol(y)")
  if ((pbirth + pdeath) > 1) stop("pbirth + pdeath must be <=1")
  if (sampler == "LIT") {
    psum= (pbirth + pdeath)
    pbirth = pbirth / psum
    pdeath = pdeath / psum
  }

  if (!is.vector(bounds_LIT)) stop("bounds_LIT must be a vector")
  if (!all(c("lbound_death", "ubound_death", "lbound_birth", "ubound_birth") %in% names(bounds_LIT))) stop("bounds_LIT must be a named vector containing 'lbound_death', 'ubound_death', 'lbound_birth', 'ubound_birth'")
  
  samplerPars= list(sampler, as.integer(niter), as.integer(burnin), as.integer(updates_per_iter), as.integer(updates_per_column), as.double(pbirth), as.double(pdeath), as.double(prob_global), as.double(tempering), as.double(truncratio), global_proposal, as.double(bounds_LIT["lbound_death"]), as.double(bounds_LIT["ubound_death"]), as.double(bounds_LIT["lbound_birth"]), as.double(bounds_LIT["ubound_birth"]), as.integer(ifelse(verbose,1,0)))
  names(samplerPars)= c('sampler','niter','burnin','updates_per_iter','updates_per_column','pbirth','pdeath','prob_global','tempering','truncratio','global_proposal','lbound_death','ubound_death','lbound_birth','ubound_birth','verbose')
  return(samplerPars)
}

#Multivariate normal density given the precision matrix
dmvnorm_prec <- function(x, sigmainv, logdet.sigmainv, mu = rep(0, ncol(sigmainv)), log = FALSE) {
  if (missing(logdet.sigmainv)) logdet.sigmainv= determinant(sigmainv,logarithm=TRUE)$modulus
  d= mahalanobis(x, center=mu, cov=sigmainv, inverted=TRUE)
  ans= -0.5 * (d + ncol(sigmainv) * log(2*pi) - logdet.sigmainv)
  if(!log) ans= exp(ans)
  return(ans)
}

#Initial estimate of precision matrix
initialEstimateGGM= function(y, Omegaini) {
    
  if (is.character(Omegaini)) {
    ans= initGGM(y, Omegaini)
  } else {
    if (ncol(Omegaini) != nrow(Omegaini)) stop("Omegaini must be a square matrix")
    if (ncol(y) != ncol(Omegaini)) stop("ncol(Omegaini) must be equal to ncol(y)")
    if (inherits(Omegaini, "matrix")) {
    ans= Matrix::Matrix(Omegaini, sparse=TRUE)
  } else if (inherits(Omegaini, c("dgCMatrix", "ddiMatrix", "dsCMatrix"))) {
    ans= Omegaini
  } else stop("Invalid Omegaini. It must be of class matrix, dgCMatrix, dsCMatrix or ddiMatrix")
  }
  return(ans)
    
}



initGGM= function(y, Omegaini) {
  
  maxiter= ifelse (ncol(y) <=100, 5, 1)  #use 1 iteration for large p, to avoid excessive init time
  
  if (Omegaini=='null') {
    
    ans= Matrix::sparseMatrix(seq_len(ncol(y)), seq_len(ncol(y)), x=rep(1,ncol(y)), dims=c(ncol(y),ncol(y)))
    
  } else if (Omegaini %in% c('glasso-bic','glasso-ebic')) {
    if (Omegaini == 'glasso-bic') {
      method <- "BIC"
      #gamma = 0
    } else if (Omegaini == 'glasso-ebic') {
      #gamma = 0.5
      method <- "EBIC"
    }
    
    sfit= huge(y, method="glasso", scr=TRUE, verbose=FALSE, nlambda = 10)
    ebic= glasso_getBIC(y=y, sfit=sfit, method=method)
    
    #topbic= which.min(sfit2$ebic.score)
    topbic= which.min(ebic)
    found= (topbic != 1) & (topbic != 10)
    
    niter= 1
    while ((!found) && (niter<=maxiter)) {
      ## huge wants a decreasing sequence
      if((niter == 1) && (topbic == 1)){## No need for lambda bigger
        break
      }
      if (topbic==10) {
        rholist= seq(sfit$lambda[10], sfit$lambda[10]/10, length=10)
      } else {
        rholist= seq(10*sfit$lambda[1], sfit$lambda[1], length=10)
      }
      
      sfit= huge(y, lambda=rholist, method="glasso", scr=TRUE, verbose=FALSE)
      ebic= glasso_getBIC(y=y, sfit=sfit, method=method)
      
      topbic= which.min(ebic)
      found= (topbic != 1) & (topbic != 10)
      niter= niter+1
      
    }
    
    ans= Matrix::Matrix(sfit$icov[[topbic]], sparse=TRUE)
    
  } else {
    stop("Omegaini must be 'null', 'glasso-ebic' or 'glasso-bic'")
  }
  
  return(ans)
  
}


#Extraction of BIC/EBIC from object returned by huge
glasso_getBIC = function(y, sfit, method='EBIC') {
  n= nrow(y); d= ncol(y)
  if (method == 'BIC') {
    gamma = 0
  } else if (method == 'EBIC') {
    gamma = 0.5
  }
  #ans = -n * sfit$loglik + log(n) * sfit$df + 4 * gamma * log(d) * sfit$df
  #Manual implementation
  logl= npar= double(length(sfit$icov))
  for (i in seq_len(length(logl))) {
    logl[i]= sum(dmvnorm_prec(y, sigmainv=sfit$icov[[i]], log=TRUE))
    npar[i]= (sum(sfit$icov[[i]] != 0) - d)/2
  }
  if (method == 'BIC') {
    ans= -2*logl + npar * log(n)
  } else if (method == 'EBIC') {
    ans= -2*logl + npar * (log(n) + 4 * gamma * log(d))
  }
  return(ans)
}




#Model log-joint (log marginal likelihood + log model prior) for the non-zero entries in colsel of Omega
#This is an R version of the C++ function GGM_rowmarg, build for debugging purposes
GGM_rowmargR= function(y, colsel, Omega, model= (Omega[,colsel]!=0), priorCoef=normalidprior(tau=1), priorModel=modelbinomprior(1/ncol(y)), priorDiag=exponentialprior(lambda=1), center=TRUE, scale=TRUE) {
  #Format input
  y = scale(y, center=center, scale=scale)
  prCoef= formatmsPriorsMarg(priorCoef=priorCoef, priorVar=priorDiag)
  prCoef= as.list(c(priorlabel=prCoef$priorCoef@priorDistr, prCoef[c('prior','tau','lambda')]))
  prModel= as.list(c(priorlabel=priorModel@priorDistr, priorPars= priorModel@priorPars))
  lambda= prCoef[["lambda"]]
  tau= prCoef[["tau"]]
  if (sum(model)>1) {
    S= t(y) %*% y
    #Obtain Omegainv, select submatrix for model
    Omegainv= solve(Omega[-colsel,-colsel,drop=FALSE])
    Omegainv_model= Omegainv[model[-colsel], model[-colsel], drop=FALSE]
    #Obtain posterior mean m and covariance Uinv
    U= (S[colsel,colsel] + lambda) * Omegainv_model + diag(1/tau, nrow=nrow(Omegainv_model))
    Uinv= solve(U)
    idx= which(model)
    idx= idx[idx != colsel]
    s= S[idx, colsel, drop=FALSE]
    m= Uinv %*% s
    #Log marginal likelihood
    logmarg= 0.5 * t(m) %*% U %*% m - 0.5 * sum(model) * log(tau) - 0.5 * log(det(U))
  } else {
    logmarg= - 0.5 * sum(model) * log(tau)
    Omegainv= Omegainv_model= Uinv= m= NULL
  }
  #Log model prior
  p= 1/ncol(y); nsel= sum(model[-colsel])
  logprior= nsel * log(p) + (length(model) - 1 - nsel) * log(1-p)
  #Returns
  ans= list(logjoint=logmarg+logprior, Omegainv=Omegainv, Omegainv_model=Omegainv_model, m=m, Uinv=Uinv)
  return(ans)
}




## huge

huge <- function (x, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, 
    method = "mb", scr = NULL, scr.num = NULL, cov.output = FALSE, 
    sym = "or", verbose = TRUE) {
    
    gcinfo(FALSE)
    est = list()
    est$method = method

    if (method == "glasso") {
        fit = huge.glasso(x, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, 
            lambda = lambda, scr = scr, cov.output = cov.output, 
            verbose = verbose)
        est$path = fit$path
        est$lambda = fit$lambda
        est$icov = fit$icov
        est$df = fit$df
        est$sparsity = fit$sparsity
        est$loglik = fit$loglik
        if (cov.output) 
            est$cov = fit$cov
        est$cov.input = fit$cov.input
        est$cov.output = fit$cov.output
        est$scr = fit$scr
        rm(fit)
        gc()
    } else {
        stop ("method not implemented")
    }

    est$data = x
    rm(x, scr, lambda, lambda.min.ratio, nlambda, cov.output, 
        verbose)
    gc()
    return(est)

}



#' The graphical lasso (glasso) using sparse matrix output 
#' 
#' @param x There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension).
#' @param lambda A sequence of decreasing positive numbers to control the regularization when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, or the thresholding in \code{method = "ct"}. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. When \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, use with care - it is better to supply a decreasing sequence values than a single (small) value.
#' @param nlambda The number of regularization/thresholding parameters. The default value is \code{30} for \code{method = "ct"} and \code{10} for \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}.
#' @param lambda.min.ratio If \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, it is the smallest value for \code{lambda}, as a fraction of the upperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "ct"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, and 0.05 \code{method = "ct"}.
#' @param scr If \code{scr = TRUE}, the lossy screening rule is applied to preselect the neighborhood before the graph estimation. The default value is  \code{FALSE}. NOT applicable when \code{method = "ct"}, {"mb"}, or {"tiger"}.
#' @param cov.output If \code{cov.output = TRUE}, the output will include a path of estimated covariance matrices. ONLY applicable when \code{method = "glasso"}. Since the estimated covariance matrices are generally not sparse, please use it with care, or it may take much memory under high-dimensional setting. The default value is \code{FALSE}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
huge.glasso = function(x, lambda = NULL, lambda.min.ratio = NULL, nlambda = NULL, scr = NULL, cov.output = FALSE, verbose = TRUE){

  gcinfo(FALSE)
  n = nrow(x)
  d = ncol(x)
  cov.input = isSymmetric(x)
  if(cov.input)
  {
    if(verbose) cat("The input is identified as the covariance matrix.\n")
    S = x
  }
  else
  {
    x = scale(x)
    S = cor(x)
  }
  rm(x)
  gc()
  if(is.null(scr)) scr = FALSE
  if(!is.null(lambda)) nlambda = length(lambda)
  if(is.null(lambda))
  {
    if(is.null(nlambda))
      nlambda = 10
    if(is.null(lambda.min.ratio))
      lambda.min.ratio = 0.1
    lambda.max = max(max(S-diag(d)),-min(S-diag(d)))
    lambda.min = lambda.min.ratio*lambda.max
    lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  }

  fit = .Call("_modelSelection_hugeglasso",S,lambda,scr,verbose,cov.output)

  fit$scr = scr
  fit$lambda = lambda
  fit$cov.input = cov.input
  fit$cov.output = cov.output

  rm(S)
  gc()
  if(verbose){
       cat("\nConducting the graphical lasso (glasso)....done.                                          \r")
       cat("\n")
      #flush.console()
  }
  return(fit)
}




## PRIOR ELICITATION ROUTINES



# Obtain prior parameters s.t. separable random matrix Theta is positive-definite with prob 0.95
# with theta[i,j] = 0 with prob 1 - edge_prob and theta[i,j] ~ N(0, sigma^2) otherwise
#
# - p : matrix dimension
# - edge_prob: probability that off-diagonal Theta[i,j] is non-zero
# - diag_distrib: if "fixed" then Theta[i,i] = diag_mean. If 'Exp' then Theta[i,i] ~ Exp(rate= 1/diag_mean)
# - diag_mean: mean of the diagonal Theta[i,i]. By default it's set to 1 for diag_distrib == "fixed", and such that P(theta_{ii} >= 1) =0.95 for diag_distrib == "Exp"
# - nsims: number of simulations
#
# Analytical lower-bound. Set sigma such that P(thetamin > 2.01 sigma sqrt(edge_prob p)) = 0.95, where thetamin is the minimum diagonal element.
# If diag_distrib == 'Exp' with theta[i,i] ~ Exp(rate= 1/diag_mean[i]), then thetamin ~ Exp(rate= sum 1 / diag_mean)
# If diag_distrib == 'fixed' then thetamin = min(diag(mean))
#
# Output
# - diag_mean: equals the input parameter if it was non-missing. Otherwise, to its default value specified above
# - offdiag_variance: sigma^2, i.e. variance of non-zero off-diagonal Theta[i,j] ~ N(0, edge_variance)
randompdmatrix_setpars <- function(p, edge_prob, diag_distrib = 'Exp', diag_mean, method='analytical', nsims=1000) {
  # Check input parameters
  if (!method %in% c('analytical','MC')) stop("method must be 'analytical' or 'MC'")
  diag_mean <- checkpars_pdmatrix_prob(p=p, diag_distrib=diag_distrib, diag_mean=diag_mean, method=method)
  # Obtain initial guess based on analytical lower-bound
  if (diag_distrib == 'Exp') {
    fexp= function(v) return((log(0.95) - pexp(2.01 * sqrt(v * edge_prob * p), rate= sum(1 / diag_mean), lower.tail=FALSE, log.p=TRUE))^2)
    #fexp= function(v) return((0.95 - pexp(2.01 * sqrt(v * edge_prob * p), rate= sum(1 / diag_mean), lower.tail=FALSE))^2)
    offdiag_variance <- optim(par=1, fn=fexp, lower=0, upper=1, method="Brent")$par
  } else {
    offdiag_variance <- min(diag_mean)^2 / (2.01^2 * p * edge_prob)
  }
  # If method != 'analytical', estimate via Monte Carlo
  if (method != 'analytical') {
    if (p > 200) {
      offdiag_variance <- randompdmatrix_setpars_extrapolate(p=p, edge_prob=edge_prob, diag_distrib=diag_distrib, diag_mean=diag_mean, nsims=nsims)$offdiag_variance
    } else {
      logvseq <- seq(log(offdiag_variance/4), log(10000 * offdiag_variance), length=nsims)  #SD from 1/2 to 100 times the analytical value
      prob_pd <- priorprob_pdmatrix(p=p, edge_prob=edge_prob, offdiag_variance=exp(logvseq), diag_mean=diag_mean, diag_distrib=diag_distrib, nsims=1)
      b <- coef(glm(prob_pd ~ logvseq, family=binomial()))
      vlow  <- -(log(1 / 0.975 - 1) + b[1])/b[2]
      vup  <- -(log(1 / 0.925 - 1) + b[1])/b[2]
      #vguess0  <- -(log(1 / 0.95 - 1) + b[1])/b[2]
      logvseq <- seq(vlow, vup, length=nsims)  #SD within range of predicted prob of positive-def in [0.927, 0.975]
      prob_pd <- priorprob_pdmatrix(p=p, edge_prob=edge_prob, offdiag_variance=exp(logvseq), diag_mean=diag_mean, diag_distrib=diag_distrib, nsims=1)
      b <- coef(glm(prob_pd ~ logvseq, family=binomial()))
      offdiag_variance  <- as.double(exp(-(log(1 / 0.95 - 1) + b[1])/b[2]))
    }
  }
  ans <- list(diag_mean=diag_mean, offdiag_variance=offdiag_variance)
  return(ans)
}

# Same as randompdmatrix_setpars, using bias-corrected analytical guess of sigma to extrapolate to large p
randompdmatrix_setpars_extrapolate <- function(p, edge_prob, diag_distrib = 'Exp', diag_mean, nsims=1000) {
  diag_mean <- min(diag_mean)
  # Obtain analytical and Monte Carlo estimates of sigma2
  pseq = round(seq(10, 80, length=10))
  sigma2_an = sigma2_mc = double(length(pseq))
  for (i in 1:length(pseq)) {
    sigma2_an[i] <- randompdmatrix_setpars(p=pseq[i], edge_prob=edge_prob, diag_mean=diag_mean, diag_distrib=diag_distrib, method="analytical")$offdiag_variance
    sigma2_mc[i] <- randompdmatrix_setpars(p=pseq[i], edge_prob=edge_prob, diag_mean=diag_mean, diag_distrib=diag_distrib, method="MC", nsims=nsims)$offdiag_variance
  }
  # Estimate bias of the analytical estimate
  b <- coef(lm(log(sigma2_an) - log(sigma2_mc) ~ pseq))
  bias_pred <- as.double(b[1] + b[2] * p)
  # Obtain bias-corrected analytical estimate
  ans <- randompdmatrix_setpars(p=p, edge_prob=edge_prob, diag_mean=diag_mean, diag_distrib=diag_distrib, method="analytical")
  ans$offdiag_variance <- exp(log(ans$offdiag_variance) - bias_pred)
  return(ans)
}


# Check input parameters of randompdmatrix_setpars, and return default diag_mean when missing
checkpars_pdmatrix_prob <- function(p, diag_distrib, diag_mean, method) {
  # Check input parameters
  if (! diag_distrib %in% c('fixed','Exp')) stop("diag_distrib muse be 'fixed' or 'Exp'")
  if (missing(diag_mean)) {
    if (diag_distrib == 'fixed') {
      diag_mean <- rep(1, p) 
    } else {
      # Since theta_{ii} ~ Exp(1/diag_mean), set diag_mean such that P(theta_{ii} >= 1)= 0.99
      diag_mean <- rep(1/0.01, p) # for 0.99 tail prob
      #diag_mean <- rep(1/0.052, p) # for 0.95 tail prob
      # Since min_i theta_{ii} ~ Exp(p/diag_mean), setting diag_mean/p= 1/0.01 achieves P(min_i theta_{ii} >= 1) = 0.99
      # diag_mean <- rep(p/0.01, p)
    }
  } else {
    if (length(diag_mean) == 1) {
      diag_mean = rep(diag_mean, p)
    } else if (length(diag_mean) != p) {
      stop("diag_mean must have length equal 1 or p")
    }
  }
  return(diag_mean)
}


# Estimate the probability that a p x p separable symmetric random matrix Theta is positive-definite
# - p : matrix dimension
# - edge_prob: probability that off-diagonal Theta[i,j] is non-zero
# - offdiag_variance: variance of non-zero off-diagonal Theta[i,j] ~ N(0, edge_variance)
# - diag_mean: mean of the diagonal Theta[i,i]
# - diag_distrib: if "fixed" then Theta[i,i] = diag_mean. If 'Exp' then Theta[i,i] ~ Exp(rate= 1/diag_mean)
# - nsims: number of simulations
# Output: Monte Carlo estimate of the probability that Theta is positive-definite
priorprob_pdmatrix <- function(p, edge_prob, offdiag_variance, diag_mean, diag_distrib = 'Exp', nsims=5000) {
  diag_mean <- checkpars_pdmatrix_prob(p=p, diag_distrib=diag_distrib, diag_mean=diag_mean)
  ans <- double(length(offdiag_variance))
  for (j in 1:length(ans)) {
    # Format input parameters
    if (length(diag_mean) == 1) {
      diag_mean = rep(diag_mean, p)
    } else if (length(diag_mean) != p) {
      stop("diag_mean must have length equal 1 or p")
    }
    if (! diag_distrib %in% c('fixed','Exp')) stop("diag_distrib muse be 'fixed' or 'Exp'")
    # Simulate Theta
    n_offdiag <- p * (p-1) / 2
    Theta <- matrix(NA, nrow=p, ncol=p)
    if (diag_distrib == 'fixed') diag(Theta) = diag_mean
    pd <- logical(nsims)
    for (i in 1:nsims) {
      # Generate diagonal entries
      if (diag_distrib == 'Exp') {
        diag(Theta) = rexp(p, rate= 1/diag_mean)
      }
      # Generate off-diagonal
      nonzero <- (runif(n_offdiag) < edge_prob)
      Theta[upper.tri(Theta)][nonzero] <- rnorm(sum(nonzero), 0, sd=sqrt(offdiag_variance[j]))
      Theta[upper.tri(Theta)][!nonzero] <- 0
      Theta[lower.tri(Theta)] <- t(Theta)[lower.tri(Theta)]
      # Check positive-definiteness
      pd[i] <- (min(eigen(Theta)$values) > 0)
    }
    ans[j] <- mean(pd)
  }
  return(ans)
}




