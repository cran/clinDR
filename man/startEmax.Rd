\name{startEmax}
\alias{startEmax}

\title{Compute starting parameter values for the 3- or 4- Emax model. }
\description{
Compute starting parameter values for iterative procedures for
estimating parameters of the 3- or 4- parameter Emax model
}
\usage{
startEmax(y,
          dose,
          baseline,
          count=rep(1,length(y)),
          modType=3,
          binary=FALSE,
          lbED50=doselev[2]/10,
          ubED50=max(doselev),
          lbLambda=0.5,
          ubLambda=5)
}
\arguments{
  \item{y}{Outcome (response) variable for the Emax modeling.}
  \item{binary}{The default is continuous (\code{binary=FALSE}).  When
  (\code{binary=TRUE}), \code{y} must be 0/1 and starting values are returned 
  for an Emax model on the logit scale.}
  \item{dose}{Dose variable corresponding to each outcome value.}
  \item{baseline}{Optional baseline covariate of same length as y.
When baseline is specified, starting values are created from anacova
adjusted dose group means.}
  \item{count}{Counts for the number of patients with each dose/y value.  
  Default is 1 (ungrouped data).}
  \item{modType}{modType=3 (default) for the 3-parameter Emax model.
  modType=4 for the 4-parameter Emax model.} 
  \item{lbED50}{If the starting ED50 is below lbED50, it is set to
lbED50.}
  \item{ubED50}{If the starting ED50 is above ubED50, it is set to
ubED50.}
  \item{lbLambda}{If the starting lambda is below lbLambda, it is set to
lbLambda.}
  \item{ubLambda}{If the starting lambda is above ubLambda, it is set to
ubLambda.}
}

\value{

Returns a vector with named elements for the starting values for a 3
or 4 parameter Emax model. The order is log(ED50), (lambda, 4 parm), 
emax, and e0. If baseline is specified, a 'beta'
starting parameter is also returned at the end of the vector.
}
\note{
The method is modified from functions created by J. Rogers and start
functions supplied with R (SSfp1).  The ED50 (and lambda) are computed using
the logit-linear relationship between the proportion of the mean
response out of the max response and the log(dose).  The method assumes
placebo data are present, but it will return a starting value even if it is
not present.  A miniumum of four dose levels is required for 4-parameter 
starting values.
}

\author{Neal Thomas }

\seealso{\code{\link{nls}},  \code{\link{emaxalt}} }
\examples{

data("examples14")
exdat<-examples14[[6]]
startEmax(exdat$y,exdat$dose)

}
\keyword{nonlinear}

