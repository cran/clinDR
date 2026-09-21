\name{targetBeta}
\alias{targetBeta}

\title{Find a scaled Beta distribution matching specified probabilities} 
\description{
Find the (a,b) parameters of a scaled Beta distribution with specified 
cummulative probabilities for two specified points from the distribution.
}
\usage{
targetBeta(minval,pminV,pmaxV,maxval=1,aInit=1,bInit=1,upB=1)
}

\arguments{
  \item{minval}{The minimum value with a targetted cummulative probability} 
  \item{pminV}{The targetted cummulative probability less than \code{minval}}
  \item{pmaxV}{The targetted cummulative probability less than \code{maxval}}
  \item{maxval}{The maximum value with a targetted cummulative probability} 
  \item{aInit}{An initial guess for the first parameter of the scaled Beta distribution with
  the specified probabilities.}
  \item{bInit}{An initial guess for the second parameter of the scaled Beta distribution with
  the specified probabilities.}
  \item{upB}{The upper limit of the scaled Beta distribution.  It is 
  specified by the user.}
}
\details{
The Beta distribution with the targetted probabilities is found from
starting values using the \code{optim} function.
}
\value{
Returns the (a,b) parameters of the scaled beta distribution if one
with the specified probabilities can be found.  An error message is 
returned otherwise.
}


\author{ Neal Thomas }
\examples{
### set quartiles at .15 and 1.0 for a beta distribution on (0,3)
targetBeta(minval=.15,pminV=0.25,pmaxV=0.75,maxval=1.0,upB=3)
}
\keyword{nonlinear}
