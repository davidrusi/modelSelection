\name{bestBIC}
\alias{bestAIC}
\alias{bestBIC}
\alias{bestEBIC}
\alias{bestIC}
\alias{bestAIC_fast}
\alias{bestBIC_fast}
\alias{bestEBIC_fast}
\alias{bestIC_fast}

\title{
Model with best AIC, BIC, EBIC or other general information criteria (getIC)
}
\description{
  Search for the regression model attaining the best value of the
  specified information criterion
}
\usage{
  bestAIC(...)

  bestAIC_fast(..., fastmethod="all")

  bestBIC(...)

  bestBIC_fast(..., fastmethod="all")

  bestEBIC(...)

  bestEBIC_fast(..., fastmethod="all")

  bestIC(..., penalty)

  bestIC_fast(..., penalty, fastmethod="all")
}
\arguments{
  \item{...}{Arguments passed on to \code{modelSelection}. The first and
    main argument is a model formula, see the examples}
  \item{penalty}{General information penalty. For example, since the AIC
    penalty is 2, bestIC(...,penalty=2) is the same as bestAIC(...)}
  \item{fastmethod}{Method used for fast model search. Set "L0Learn" to
    use the L0Learn package, "L1" to use LASSO (glmnet package), 
    "adaptiveL1" for adaptive LASSO, and
    "CDA" for coordinate descent. Not all these options may be
    available for some GLM families (e.g. L0Learn only supports
    \code{family="normal"} and \code{"binomial"}). The default "all"
    uses all available methods}
}
\value{Object of class \code{icfit}. Use (coef, summary,
  confint, predict) to get inference for the top model,
  and \code{help(icfit-class)} for more details on the returned object.
}
\details{

  \code{bestAIC}, \code{bestBIC}, \code{bestEBIC}, \code{bestIC} perform
  full model enumeration when possible and otherwise resort to MCMC to
  explore the models, as discussed in function \code{modelSelection}.
  
  \code{bestAIC_fast}, \code{bestBIC_fast}, \code{bestEBIC_fast},
  \code{bestIC_fast} use a faster algorithm. It first identifies
  a subset of promising models, and then computes the specified
  criterion for each of them to find the best one within the subset.
  For Gaussian and binary outcomes it uses function \code{L0Learn.fit} from
  package \code{L0Learn} (Hazimeh et al, 2023), which combines
  coordinate descent with local combinatorial search to find good models
  of each size.
  L1 returns all the models found in the LASSO regularization path.
  CDA returns a single model found by coordinate descent,
  i.e. adding/dropping one covariate at a time to improve the
  specified criterion (BIC, AIC, ...). 
  
  \code{bestBIC} and the other functions documented here take similar
  arguments to those of \code{modelSelection}, but here no priors
  on models or parameters are needed.
  
  Let p be the total number of parameters and n the sample size. The BIC of a model k with p_k
  parameters is

  - 2 L_k + p_k log(n)

  the AIC is
  
  - 2 L_k + p_k 2
  
  the EBIC is

  - 2 L_k + p_k log(n) + 2 log(p choose p_k)

  and a general information criterion with a given model size penalty

  - 2 L_k + p_k penalty


  The MCMC model search is based on assigning a probability to each
  model, and then using MCMC to sample models from this
  distribution. The probability of model k is

  exp(- IC_k / 2) / sum_l exp(- IC_l / 2)

  where IC_k is the value of the information criterion (BIC, EBIC...)

  Hence the model with best (lowest) IC_k has highest probability, which
  means that it is likely to be sampled by the MCMC algorithm.
 
}
\author{
David Rossell
}
\seealso{
\code{\link{modelSelection}} to perform Bayesian model selection, and
for a description of the arguments that should be passed in \code{...}
}
\references{
H. Hazimeh, R. Mazumder, T. Nonet. L0learn: A scalable package for sparse learning using l0 regularization. Journal of Machine Learning Research 24.205 (2023): 1-8.
}
\examples{
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x \%*\% theta + rnorm(100)
ybin <- y>0
df <- data.frame(y, ybin, x)

#BIC for all models (the intercept is also selected in/out)
fit= bestBIC(y ~ X1 + X2, data=df)
fit

#Same, but setting the BIC's log(n) penalty manually
#change the penalty for other General Info Criteria
#n= nrow(x)
#fit= bestIC(y ~ X1 + X2, data=df, penalty=log(n))

summary(fit) #usual GLM summary

coef(fit) #MLE under top model

#confint(fit) #conf int under top model (requires MASS package)


#Binary outcome
fit2= bestBIC(ybin ~ X1 + X2, data=df, family='binomial')
fit2


}
\keyword{ models }
\keyword{ distribution }
