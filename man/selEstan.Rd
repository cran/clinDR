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

selEstan(emod=c('basemodel.rds','mrmodel.rds'))
}

\arguments{
\item{emod}{Two parameterizations of the emax function are currently supported.  
'basemodel'  uses the maximal effect 'emax' parameter.  'mrmodel'  uses the effect
of the drug at a high dose specified by the user versus placebo.  The 'emax' effect model is
deprecated and will be eliminated.}
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
	estan<-selEstan()
}
}

