##############################################################################################
##
## ROUTINES TO COMPUTE INTEGRATED LIKELIHOODS
##
##############################################################################################

marginalLikelihood <- function(
  sel, y, x, data, smoothterms, nknots=9, groups=1:ncol(x), family="normal",
  priorCoef, priorGroup, priorVar=igprior(alpha=0.01,lambda=0.01),
  priorSkew=momprior(tau=0.348), neighbours, phi, method='auto', adj.overdisp='intercept', hess='asymp',
  optimMethod, optim_maxit, initpar='none', B=10^5, logscale=TRUE, XtX, ytX
) {
  #Check input
  # format input data
  tmp <- formatInputdata(y=y,x=x,data=data,smoothterms=smoothterms,nknots=nknots,family=family)
  x <- tmp$x; y <- tmp$y; is_formula <- tmp$is_formula
  splineDegree <- tmp$splineDegree
  if (!is.null(tmp$groups)) groups <- tmp$groups
  hasgroups <- tmp$hasgroups
  if (!is.null(tmp$constraints)) constraints <- tmp$constraints
  outcometype <- tmp$outcometype; uncens <- tmp$uncens; ordery <- tmp$ordery
  typeofvar <- tmp$typeofvar
  p= ncol(x); n= length(y)
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  sumy2 <- as.double(sum(y^2)); sumy <- as.double(sum(y))
  colsumsx <- as.double(colSums(x))
  #
  familyint= formatFamily(family, issurvival= length(uncens)>0)$familyint
  if (familyint == 22) { sumlogyfact= as.double(sum(lgamma(y+1))) } else { sumlogyfact= as.double(0) } #Poisson regression
  # check prior and set defaults if necessary
  if (missing(priorCoef)) {
      defaultprior= defaultmom(outcometype=outcometype,family=family)
      priorCoef= defaultprior$priorCoef; priorVar= defaultprior$priorVar
  }
  if (missing(priorGroup)) { if (length(groups)==length(unique(groups))) { priorGroup= priorCoef } else { priorGroup= groupzellnerprior(tau=n) } }

  if (missing(phi)) { knownphi <- as.integer(0); phi <- double(0) } else { knownphi <- as.integer(1); phi <- as.double(phi) }

  if (!missing(neighbours)) { Dmat= icar_dmatrix(neighbours) } else { Dmat= diag(p) }

  # format arguments for .Call
  thinit= getthinit(y=y, x=x, family=family, initpar=initpar, enumerate=TRUE)
  usethinit= thinit$usethinit; thinit= thinit$thinit

  method <- formatmsMethod(method=method, usethinit=usethinit, optimMethod=optimMethod, optim_maxit=optim_maxit, priorCoef=priorCoef, priorGroup=priorGroup, knownphi=0, outcometype=outcometype, family=family, hasgroups=hasgroups, adj.overdisp=adj.overdisp, hess=hess)
  optimMethod <- method$optimMethod; optim_maxit <- method$optim_maxit; adj.overdisp <- method$adj.overdisp; hesstype <- method$hesstype; method <- method$method
  #hesstype <- as.integer(ifelse(hess=='asympDiagAdj',2,1)); optimMethod <- as.integer(ifelse(optimMethod=='CDA',2,1))
    
  B <- as.integer(B)
  tmp= codeGroupsAndConstraints(p=p,groups=groups)
  ngroups= tmp$ngroups; constraints= tmp$constraints; invconstraints= tmp$invconstraints; nvaringroup=tmp$nvaringroup; groups=tmp$groups
  tmp= formatmsPriorsMarg(priorCoef=priorCoef, priorGroup=priorGroup, priorVar=priorVar, priorSkew=priorSkew, n=n)
  r= tmp$r; prior= tmp$prior; priorgr= tmp$priorgr; tau=tmp$tau; taugroup=tmp$taugroup; alpha=tmp$alpha; lambda=tmp$lambda; taualpha=tmp$taualpha; fixatanhalpha=tmp$fixatanhalpha
  a= tmp$a

  if (!is_formula) {
    sel <- check_sel_groups(sel, groups)
    sel <- as.integer(sel-1); nsel <- as.integer(length(sel))
  } else {
    if (!missing(sel)) warning("y is of type formula: ignoring sel argument")
    sel <- as.integer(seq(p)-1)
    nsel <- length(sel)
  }

  ans <- marginalLikelihoodCI(knownphi, sel, nsel, familyint, prior, priorgr, n, p, y, uncens, sumy2, sumy, sumlogyfact, x, colsumsx, XtX, ytX, method, adj.overdisp, hesstype, optimMethod, optim_maxit, thinit, usethinit, B, alpha, lambda, tau, taugroup, taualpha, fixatanhalpha, r, a, groups, ngroups, nvaringroup, constraints, invconstraints, Dmat, logscale)
  return(ans)
}

check_sel_groups <- function(sel, groups) {
  p <- length(groups); seqp <- seq(p)
  if (is.logical(sel)) sel <- seqp[sel]
  if (any(sel > p)) stop("found index in sel larger than ncol(x). Please make sure all indexes refer to existing variables")
  invsel <- seqp[!(seqp %in% sel)]
  if (any(groups[sel] %in% groups[invsel])) stop("selected indexes incompatible with defined groups. Make sure each group is selected or discarded at once")
  return(sel)
}


##############################################################################################
## EXPECTATION OF PRODUCT OF SQUARES
##############################################################################################


eprod <- function(m, S, power=1, dof= -1) {
  #Mean of prod (x_i)^(2*power) when x_i ~ T_dof(mu,sigma). Set dof=-1 for N(mu,sigma). Written by John Cook
  ans <- eprod_I(as.double(m),as.double(S), as.integer(length(m)), as.integer(power), as.double(dof))
  ans
}

