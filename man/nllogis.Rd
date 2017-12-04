\name{nllogis}
\alias{nllogis}
\title{The negative log likelihood function for a 3- or 4-
parameter Emax model on the logit scale for binary dose response.
}
\description{
The negative log likelihood function evaluated with a single input set of
parameters for the binary Emax model on the logistic scale.  For use
with function \code{\link{fitEmax}}
}
\usage{
nllogis(parms,y,dose,
          prot=rep(1,length(y)),
          count=rep(1,length(y)))
}

\arguments{
  \item{parms}{Emax model parameter values.  The
order of the variables is (log(ED50),Emax,E0) or
(log(ED50),lambda,Emax,E0). There must be an E0 for each protocol.  
Note the transformation of ED50.}
  \item{y}{Binary outcome variable for each patient.  Missing values
are deleted. Must be coded 0/1.}
  \item{dose}{Dose for each patient}
  \item{prot}{Protocol (group) membership used to create multiple intercepts.  
  The default is a single protocol.  The value of \code{prot} must be 1,2,3,..}
  \item{count}{Counts for the number of patients with each dose/y value.  
  Default is 1 (ungrouped data).}
}

\details{
The negative log likelihood for the 3- or 4- Emax
model on the logit scale for binary data.  Note the ordering of the parameters
and their transformations.  A 3 vs 4 parameter model is deterimined by
the length of \code{parms}.}

\value{
Negative log likelihood value is returned.
}


\author{Neal Thomas }

\seealso{\code{\link{nlm}},  \code{\link{fitEmax}} }

\examples{
data(examples14)
with(examples14[[8]],nllogis(parms=c(log(.17),-3.26,-0.15), y, dose))
}

\keyword{nonlinear}