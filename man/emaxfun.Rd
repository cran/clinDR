\name{emaxfun}
\alias{emaxfun}

\title{ Vectorized versions of the hyperbolic and sigmoidal Emax models }
\description{
Evaluate Emax models for a vector of dose levels for multiple sets
of parameters.
}
\usage{
emaxfun(dose, parm)
}

\arguments{
  \item{dose}{ A vector (or scalar) of dose levels }
  \item{parm}{A vector or matrix with columns containing
              \code{log(ed50)}, Hill parameter if sigmoid model,
              \code{emax},\code{e0}}
}

\details{
The Hill parameter is omitted from \code{parm} for the hyperbolic model
}

\value{
Returns a matrix of Emax function evaluations.  The rows correspond to the
parameter replications, and the columns correspond to the dose levels.}

\note{The ordering of the parameters was selected to faciliate use of
the 'plinear' algorithm in function nls.
}

\author{ Neal Thomas }

\seealso{\code{\link{dlogis}} }
\examples{

doselev<-c(0,5,25,50,100)
e0<-2.465375 
ed50<-67.481113 

dtarget<-100
diftarget<-9.032497
lambda=2
emax<-solveEmax(diftarget,dtarget,log(ed50),lambda,e0)

parm<-c(log(ed50),lambda,emax,e0)
plot(doselev,emaxfun(doselev,parm))

}
\keyword{nonlinear}
