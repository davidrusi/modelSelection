\name{postProb}
\alias{postProb}
\alias{postProb-methods}
\alias{postProb,msfit-method}
\alias{postProb,msfit_ggm-method}
\alias{postProb,cilfit-method}
\alias{postProb,localtest-method}
\alias{postProb,mixturebf-method}
\alias{postProbPlot}
\alias{postProbPlot-methods}
\alias{postProbPlot,msfit-method}


\title{
Obtain and plot posterior model probabilities
}
\description{
  Obtain posterior model probabilities after running Bayesian model selection
}
\usage{
  postProb(object, nmax, method='norm')

  postProbPlot(object, nmax, method='norm', sort_variables=TRUE)
}
\arguments{
  \item{object}{Object of class msfit returned by \code{modelSelection},
    class mixturebf returned by \code{bfnormmix},
    class cilfit returned by \code{cil}
    or class localtest returned by \code{localnulltest}}

  \item{nmax}{Maximum number of models to report (defaults to no max)}

  \item{method}{Only when \code{class(object)} is msfit.
    For 'norm' probabilities are obtained by renormalizing the
  stored integrated likelihoods, for 'exact' they are given by the proportion
  of MCMC visits to each model. 'norm' has less variability but can be biased
  if the chain has not converged.}

  \item{sort_variables}{If set to \code{TRUE} variables in the y-axis are 
plotted according to their posterior inclusion probabilities, 
otherwise the original order in the specified design matrix is used}

  \item{max_variables}{Maximum number of variables to be displayed in the
y-axis. If missing, all variables are displayed}

}
\value{
  A \code{data.frame} with posterior model probabilities in column pp.
  Column modelid indicates the indexes of the selected covariates (empty
  for the null model with no covariates)

  For \code{msfit_ggm} objects, a list with posterior probabilities
  (element "pp") and the corresponding model (element "modelid")
  indicating what edges are non-zero

  For \code{mixturebf} objects, posterior probabilities for the
  specified number of components

  For \code{localtest} objects, posterior probabilities of a local covariate
  effect at various regions
}
\author{
David Rossell
}
\seealso{
\code{\link{modelSelection}} to perform model selection
}
\examples{
#See help(modelSelection)
}
\keyword{ models }
\keyword{ distribution }
