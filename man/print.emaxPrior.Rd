\name{printemaxPrior}
\alias{print.emaxPrior}
\title{
Print protocol or sap text describing the prior distribution for
the model parameters of the input emaxPrior object
}
\description{
Print templated description of the prior distribution for the 
Emax model parameters.  The level of detail is adjusted for protocol/sap. 
By default, the prior object is printed as a list without documentation.
}
\usage{
\method{print}{emaxPrior}(x, doc=FALSE, diffuse=NULL, file="", 
		modType=c('4','3'), docType=c('sap','protocol'), ...)
}
\arguments{
  \item{x}{Object created by function \code{emaxPrior.control}
}
	\item{doc}{When \code{TRUE}, documentation for the prior distribution
	is returned.  Default is \code{FALSE} and the prior input is returned
	as a list for default printing}
	\item{diffuse}{When \code{TRUE}, the scale parameters are described
	as creating diffuse prior distributions for the corresponding 
	efficacy parameters.  When \code{FALSE}, sections are identified
	where the user must add justification for the informative 
	prior distributions.  An error message is printed if 
	\code{doc=TRUE} and \code{diffuse} is not specified
}
	\item{file}{File for ascii output}
	\item{modType}{Character value ('4' or '3') that determines whether
	the 4-parameter sigmoidal Emax paramter is included, or the 3-parameter
	hyperbolic model is assumed with the Hill (slope) parameter set to 1
}
	\item{docType}{When \code{'protocol'}, the prior description
	is less detailed.
}
	\item{...}{No other inputs are supported}
}
\details{
If the object is entered at the command line, the implied print
function is called without the required \code{diffuse} flag.  The
object will be printed as a list.  The list output will be
followed by error/warning messages noting the absence of the
required input.
}
\value{
Ascii text or text file that can be edited for inclusion in
a protocol/sap
}
\author{
Neal Thomas
}
\seealso{
\code{\link{emaxPrior.control}}
}
\examples{
prior<-emaxPrior.control(epmu=0,epsca=4,difTargetmu=0,
		difTargetsca=4,dTarget=20,
		p50=(2+5)/2,
		sigmalow=0.01,sigmaup=3)
				
print(prior,doc=TRUE,diffuse=TRUE)

}
\keyword{ prior }
\keyword{ Emax }

