#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()
#Rcpp::compileAttributes("~/github/mombf")

library(modelSelection)
library(mvtnorm)
set.seed(1)

## GGM examples
library(modelSelection)
library(mvtnorm)
set.seed(1)
p <- 5
n <- 5000
set.seed(1)
Omega= diag(p)
Omega[abs(col(Omega) - row(Omega))==1]= 0.5
y = rmvnorm(n = n, sigma = solve(Omega))
#
niter= 200
burnin= 0
updates_per_iter= p; updates_per_column= p
pbirth= 0.75; pdeath= 1-pbirth

priorModel= modelbbprior(alpha.p=1, beta.p=1)
fit.bb= modelSelectionGGM(y, sampler='Gibbs', niter=niter, Omegaini=Omega, priorModel=priorModel, global_proposal="none", updates_per_iter=updates_per_iter, updates_per_column=updates_per_column, pbirth= pbirth, pdeath=pdeath)
