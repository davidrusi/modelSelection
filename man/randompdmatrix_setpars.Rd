\name{randompdmatrix_setpars}
\alias{randompdmatrix_setpars}
\title{

Set parameters such that symmetric random matrix 
is positive-definite with high probability

}
\description{

Obtain prior parameters s.t. separable symmetric random matrix Theta 
is positive-definite with prob 0.95

Its off-diagonal entries are drawn as theta[i,j] = 0 with prob 1 - edge_prob 
and theta[i,j] ~ N(0, sigma^2) otherwise, independently across i < j

Its diagonal entries theta[i,i] can either be fixed to diag_mean or randomly 


}
\usage{

randompdmatrix_setpars(p, edge_prob, diag_distrib = 'Exp', diag_mean, method='analytical', nsims=1000)

}
\arguments{

\item{p}{Matrix dimension}

\item{edge_prob}{Probability that off-diagonal Theta[i,j] is non-zero}

\item{diag_distrib}{If "fixed" then Theta[i,i] = diag_mean. If "Exp" then Theta[i,i] ~ Exp(rate= 1/diag_mean)}

\item{diag_mean}{Mean of the diagonal Theta[i,i]. By default it's set to 1 for diag_distrib == "fixed", and such that P(theta_{ii} >= 1) =0.95 for diag_distrib == "Exp"}

\item{nsims}{Number of simulations}

}
\details{

The analytical estimate of the off-diagonal variance sigma 
is obtained as follows. 
If diag(Theta) is fixed, then sigma is set such that

thetamin > 2.01 sigma sqrt(edge_prob p)

where thetamin = min(diag(Theta)), since as p -> infty this guarantees
that Theta is positive-definite with probability 1 (Carter & Rossell 2026)

If theta[i,i] ~ Exp(1/diag_mean[i]), then sigma is set such that

P(thetamin > 2.01 sigma sqrt(edge_prob p)) = 0.95, 

using that thetamin ~ Exp(rate= sum 1 / diag_mean).
If the probability above goes to 1 as p -> infty, this guarantees
that Theta is positive-definite with probability 1 (Carter & Rossell 2026)

The Monte Carlo estimate of sigma is based on an iterative logistic
regression approach. Specifically, it simulates nsim matrices
Theta for a range of sigma values around its analytical estimate,
estimates P(Theta is pos-def) as a function of sigma via 
logistic regression and updates sigma such that the predicted
P(Theta is pos-def | sigma)= 0.95.
The process is then repeated: given the updated sigma, 
nsim new matrices are drawn around that sigma
and a new logistic regression is fit to update sigma

}
\value{

A list with two entries

\item{diag_mean}{Equals the input parameter if it was non-missing. 
Otherwise it's set to its default value, specified above}

\item{offdiag_variance}{sigma2, i.e. variance of non-zero 
off-diagonal Theta[i,j] ~ N(0, edge_variance)}

}
\references{
Jack Carter, David Rossell. Positive-definiteness in separable priors: effects on prior interpretability and inference. arXiV 2026
}
\author{
David Rossell
}
\seealso{
\code{
\link{modelSelectionGGM}} to perform model selection for graphical models
}
\examples{

# Prior parameters based on asymptotics as p -> infty
randompdmatrix_setpars(p=10, edge_prob=0.1, diag_distrib='Exp', method='analytical')

# Prior parameters based on Monte Carlo for fixed p
#randompdmatrix_setpars(p=10, edge_prob=0.1, diag_distrib='Exp', method='MC')

}
\keyword{ models }
\keyword{ distribution }
