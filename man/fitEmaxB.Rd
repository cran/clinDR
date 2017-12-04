\name{fitEmaxB}
\alias{fitEmaxB}
\title{
Bayesian fit of hyperbolic or sigmoidal Emax models to continuous/binary dose
response data.
}
\description{
Uses Rpackage \code{\link{rstan}} to fit a Bayesian hyperbolic or sigmoidal 
Emax model.  Different intercepts for multiple protocol-data are supported.  For binary
data, the Emax model is on the logit scale.
}
\usage{
fitEmaxB(y, dose, prior, modType = 3, prot = rep(1, length(y)), 
	count = rep(1, length(y)), binary = FALSE, msSat = NULL,
	pboAdj = FALSE, mcmc = mcmc.control(),estan = NULL, diagnostics = TRUE)
}

\arguments{

  \item{y}{Outcome for each patient.  Missing \code{Y} values
are are not permitted.  Dose/protocol group means for grouped continuous data.  For binary
data, \code{y} must be 0/1 and counts must be supplied for each 0/1 value.}
  \item{dose}{Dose for each patient.}
  \item{prior}{Prior specification through a list of parameter values.  See \code{\link{prior.control}}
  for details.
}
  \item{modType}{modType=3 (default) for the 3-parameter hyperbolic Emax model.
  modType=4 for the 4-parameter sigmoidal Emax model.} 
  \item{prot}{Protocol (group) membership used to create multiple intercepts.  
  The default is a single protocol.  The  prior disribution for the placebo
  response is re-used independently for each intercept.}
  \item{count}{Counts for the number of patients when the \code{Y} are dose 
  continuous group means or binary 0/1 values.  
  Default is 1 (ungrouped data).}
  \item{binary}{When \code{TRUE}, the \code{y} are assumed to be coded \code{0/1}, and 
  the the means reported are proportions.  The Emax model is specified on the logit 
  scale, and proportions are estimated from the model by back-transformation.}
  \item{msSat}{If continuous \code{Y} are dose/protocol group means rather than
  individual measurements, the within group variance, \code{msSat}, should be 
  supplied.  This variance is the mean square from the model saturated in
  dose and protocol.  It is used for goodness-of-fit (GOF) testing, and to 
  improve the residual variance estimate for the Emax model.  If it is not
  supplied, statistics needed for GOF will not be available, and the
  residual SD (and associated SE) will have low degrees of freedom.}
  \item{pboAdj}{For published data with only pbo-adjusted dose group means and
  SEs, the model is fit without an intercept(s).  If initial parameters
  are supplied, the intercept (E0) should be assigned \code{0}.  
  A zero for the placebo mean should not be included in \code{Y}.  
  This option is not available for binary data.  Potential correlation between between
  placebo-adjusted means is ignored.} 
  \item{mcmc}{
  Inputs controlling \code{rstan} execution.  See \code{\link{mcmc.control}}
  for details.}
  \item{estan}{The compiled \code{rstan} Emax model is usually loaded
  automatically.  It can be load to an object using the function
  \code{selEstan} and passed to \code{fitEmaxB} for repeated executions to
  improve efficiency and stability.}
  \item{diagnostics}{Print output from rstan.}
}
\details{
The function \code{compileStanModels} must be executed once to create compiled STAN
code before \code{fitEmaxB} can be used.

MCMC fit of a Bayesian hyperbolic or sigmoidal Emax model.  The prior
distributions available are based on the publication  
Thomas, Sweeney, and Somayaji (2014) and
Thomas and Roy (2016).}
\value{
A list assigned class "fitEmaxB" with:
\item{estanfit}{The \code{rstan} object with the model fit.}
\item{y, dose, prot, count, modType, binary, pboAdj, 
msSat, prior, mcmc}{Input values.}
}
\references{
Thomas, N., Sweeney, K., and Somayaji, V. (2014).  Meta-analysis of clinical dose response in a large drug development portfolio, Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317.
<doi:10.1080/19466315.2014.924876>

Thomas, N., and Roy, D. (2016).  Analysis of clinical dose-response in small-molecule drug development:  2009-2014. Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317
<doi:10.1080/19466315.2016.1256229>
}
\author{
Neal Thomas
}


\seealso{
\code{\link{fitEmax}}, \code{\link{predict.fitEmaxB}}, \code{\link{plot.fitEmaxB}}
}
\examples{
\dontrun{

data("examples14")
exdat<-examples14[[1]]

prior<-prior.control(epmu=0,epsd=4,emaxmu=0,emaxsd=4,p50=0.1,
				sigmalow=0.01,sigmaup=3)
										
mcmc<-mcmc.control(chains=3)

fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
				count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)
}
}

\keyword{ dose response}
\keyword{ mcmc }
