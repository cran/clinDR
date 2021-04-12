\name{runShiny}
\alias{runShiny}
\title{
Shiny app for function \code{emaxsim(B)} 
}
\description{Shiny app for function \code{emaxsim(B)}}
\usage{
runShiny()
}
\author{
Neal Thomas, Mike K. Smith
}
\note{
The code section of the shiny app provides the code required
for batch execution of the current shiny results.

The 'Analysis' section of the shiny app must be visited before
an example can be run.  

For Bayesian output, the \code{clinDR}
package function \code{compileStanModels()} must be executed
once before using the shiny app or any of the package functions
utilizing Bayes methods.
}
\examples{
	if (interactive()) {
	  runShiny ()
	}
}
\seealso{
\code{\link{emaxsimB}}
}


