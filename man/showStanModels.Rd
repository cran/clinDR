\name{showStanModels}
\alias{showStanModels}
\title{
Display \code{STAN} model code.
}
\description{
Display the \code{STAN} Bayesian model code for fitting Emax models
}
\usage{
showStanModels(emod=c('basemodel.stan','mrmodel.stan'))
}
\arguments{
\item{emod}{Two parameterizations of the emax function are currently supported.  
'basemodel'  uses the maximal effect 'emax' parameter.  'mrmodel'  uses the effect
of the drug at a high dose specified by the user versus placebo.  The 'emax' effect model is
deprecated and will be eliminated.}
}
\author{
Neal Thomas
}

\seealso{
\code{\link{fitEmaxB}}, \code{\link{emaxsimB}}
}

\examples{
\dontrun{
showStanModels()
}
}

