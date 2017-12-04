\name{selEstan}
\alias{selEstan}
\title{
Select a pre-compiled \code{rstan} Emax model
}
\description{
Emax models for use in \code{fitEmaxB} and \code{emaxsimB} which have been pre-compiled
are loaded for use outside of the the fiting functions.  This is most useful for
repeated simulations in which the loading of the compiled models from
a disk file can be performed once.  \code{fitEmaxB} will load the model
automatically for  single execution, so the model does not
need to be pre-loaded.
}
\usage{
selEstan(modType = 3, binary = FALSE, pboAdj = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modType}{
3- or 4- hyperbolic or sigmoidal Emax model.
}
  \item{binary}{
An Emax model on the logistic scale for binary outcomes.
}
  \item{pboAdj}{
Model omitting the intercept when only placebo-corrected data are 
available.  Only valid with continuous outcomes.
}
}

\value{
An Emax 'stanmodel'.
}

\author{
Neal Thomas
}

\seealso{
\code{\link{fitEmaxB}}, \code{\link{emaxsimB}}
}
\examples{
\dontrun{
	estan<-selEstan(modType=4)
}
}

