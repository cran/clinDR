\name{clinDR-package}
\alias{clinDR-package}
\docType{package}
\title{
Bayesian and maximum likelihood Emax model fitting, graphics and 
simulation for clinical dose response.
}
\description{
The functions \code{\link{fitEmax}} and \code{\link{fitEmaxB}} fit an Emax model
to binary or continuous data using maximum likelihood or Bayesian estimation. 
They have several generic supporting functions.
Functions to produce plots associated with dose response
analyses are (\code{\link{plotD}}, \code{\link{plotB}},
\code{\link{plot.fitEmax}},\code{\link{plot.fitEmaxB}}).   
The functions \code{\link{emaxsim}} and \code{\link{emaxsimB}} perform
simulations of 4- and 3-parameter Emax ML or Bayesian estimation. 
The ML estimates are replaced with alternative
model fits when the primary estimation fails.  Several supporting
functions are supplied to analyze the output of \code{\link{emaxsim}}
and \code{\link{emaxsimB}},
including analyses for specific simulated data sets.  All of the data sets
from dose response meta analyses are included in 
\code{\link{metaData}}.  
}
\details{
The function \code{compileStanModels} must be executed once after the package is
installed to create compiled \code{STAN} Emax models before the Bayes functions 
in the package can be executed.  This requires 3-10 minutes to complete on most
machines. The compiled code is 32-bit or 64-bit specific, and both must be created if
both versions of R are used.

The Bayesian computations use the R package \code{rstan}.  It can be installed from CRAN. Windows users should check the instructions for \code{rstan} at the \url{http://mc-stan.org} and \url{https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}. 
Note that \code{Rtools} must be installed, which is a simple, but often overlooked step. Instructions for its installation are given in the second URL.
}

\author{
Neal Thomas [aut, cre], Jing Wu[aut]
}

\keyword{ package }

\seealso{
\code{\link{DoseFinding}}} 

