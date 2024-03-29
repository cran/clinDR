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
fitEmaxB(y, dose, prior, modType = 4, prot = rep(1, length(y)), 
	count = rep(1, length(y)), xbase=NULL,
	binary = FALSE, msSat = NULL,vcest=NULL,
	pboAdj = FALSE, mcmc = mcmc.control(), estan = NULL, 
	diagnostics = TRUE, nproc = getOption("mc.cores", 1L))
}

\arguments{

  \item{y}{Outcome for each patient.  Missing \code{Y} values
are are not permitted.  Dose/protocol group means for grouped continuous data.  For binary
data, \code{y} must be 0/1 and counts must be supplied for each 0/1 value.}
  \item{dose}{Dose for each patient.}
  \item{prior}{Prior specification through an object of type 'emaxPrior' or 'prior'.  
  See \code{\link{emaxPrior.control}} and \code{\link{prior.control}} for details.  
  The 'emaxPrior' specifies
  the magnitude of the potential effect for a specified dose (typically the highest
  anticipated dose and/or the dose in a POC study), while the 'prior' specifies the
  theoretical maximum effect (the emax parameter).  The 'prior' specification is 
  deprecated and will be removed.}
  \item{modType}{modType=3 (default) for the 3-parameter hyperbolic Emax model.
  modType=4 for the 4-parameter sigmoidal Emax model.} 
  \item{prot}{Protocol (group) membership used to create multiple intercepts.  
  The default is a single protocol.  The  prior disribution for the placebo
  response is re-used independently for each intercept.}
  \item{count}{Counts for the number of patients when the \code{Y} are dose 
  continuous group means or binary 0/1 values.  
  Default is 1 (ungrouped data).}
  \item{xbase}{A matrix of baseline covariates with rows corresponding to 
  \code{y} that enter as linear additive predictors.  The baseline covariates must 
  be centered about
  their (protocol-specific) means.  \code{xbase} does not include an intercept or
  protocol indicators.  Covariates cannnot be specified with PBO adjusted or 
  aggregated input.}
  \item{binary}{When \code{TRUE}, the \code{y} are assumed to be coded \code{0/1}, and 
  the means reported are proportions.  The Emax model is specified on the logit 
  scale, and proportions are estimated from the model by back-transformation.}
  \item{msSat}{If continuous \code{Y} are dose/protocol group means rather than
  individual measurements, the within group variance, \code{msSat}, should be 
  supplied.  This variance is the mean square from the model saturated in
  dose and protocol.  It is used to 
  improve the residual variance estimate for the Emax model.  If it is not
  supplied, the
  residual SD (and associated SE) will have low degrees of freedom.}
  \item{vcest}{The input, \code{Y}, can be estimates of dose group
  responses from a first-stage
  model.  The \code{vcest} is the variance-covariance matrix of the
  model-based estimates.  The most common usage is when a saturated model
  is fit using maximum likelihood estimation to longitudinal data 
  to produce dose group estimates that are valid under the MAR assumption
  for missing values.  Other applications are possible, see the Pinheiro, 
  et al reference.  The \code{count} and \code{msSat} are ignored when
  \code{vcest} is specified.  Covariates \code{xbase} cannot be specified
  with \code{vcest}, but covariates can be included in the first stage 
  modeling.
  }
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
  \item{diagnostics}{Printed output from rstan.  See \code{Details} for more information.}
  \item{nproc}{The number of processor requested for \code{STAN} MCMC
  computations.  Defaults to the value set by the 
  \code{rstan} installation.   When set explicitly, \code{nproc} 
  is usually 1 or the number of MCMC chains. If greater than the
  number of chains, it is set to the number of chains.}
}
\details{
The function \code{compileStanModels} must be executed once to create compiled STAN
code before \code{fitEmaxB} can be used.

MCMC fit of a Bayesian hyperbolic or sigmoidal Emax model.  The prior
distributions available are based on the publication  
Thomas, Sweeney, and Somayaji (2014), 
Thomas and Roy (2016), and Wu, et al (2017).

The posterior distributions are complex because the distributions of the \code{Emax} and \code{ED50} parameters 
change substantially as a function of the \code{lambda}, often
creating 'funnel' type conditions. Small numbers of 
divergences are common with the 4-parameter model and do not 
appear easily avoided.
Extensive simulation using evaluations with \code{emaxsimB}
support the utility of the resulting approximate 
posterior distributions. The number of divergences can be
viewed using \code{diagnostics=TRUE}. The usual convergence
diagnostics should always be checked.
}
\value{
A list assigned class "fitEmaxB" with:
\item{estanfit}{The \code{rstan} object with the model fit.}
\item{
y, dose, prot,count,
nbase[rows of xbase],xbase,dimFit[rows of vcest],
vcest,modType,binary,pboAdj,msSat,
prior,mcmc,localParm[TRUE when
\code{emaxPrior} specified]		
}{Input values.}
}
\note{The default modType was changed from 3 to 4 for clinDR version >2.0}
\references{
Thomas, N., Sweeney, K., and Somayaji, V. (2014).  Meta-analysis of clinical dose response in a large drug development portfolio, Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317.
<doi:10.1080/19466315.2014.924876>

Pinheiro, J., Bornkamp, B., Glimm, E., and Bretz, F. (2014).  
Model-based dose finding under model uncertainty using general 
parametric models, Vol. 33, No. 10, 1646-1661 <doi:/10.1002/sim.6052>

Thomas, N., and Roy, D. (2016).  Analysis of clinical dose-response in small-molecule drug development:  2009-2014. Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317
<doi:10.1080/19466315.2016.1256229>

Wu, J., Banerjee, A., Jin, B. Menon, M. S., Martin, S. and Heatherington, A. (2017). Clinical dose response for a broad set of biological products: A model-based meta-analysis. Statistical Methods in Medical Research. <doi:10.1177/0962280216684528>
}
\author{
Neal Thomas
}


\seealso{
\code{\link{fitEmax}}, \code{\link{predict.fitEmaxB}}, \code{\link{plot.fitEmaxB}}, \code{\link{coef.fitEmaxB}}
}
\examples{
\dontrun{

data("metaData")
exdat<-metaData[metaData$taid==1,]

prior<-emaxPrior.control(epmu=0,epsca=4,difTargetmu=0,difTargetsca=4,dTarget=20,
				p50=(2+5)/2,
				sigmalow=0.01,sigmaup=3)
										
mcmc<-mcmc.control(chains=3)

msSat<-sum((exdat$sampsize-1)*(exdat$sd)^2)/(sum(exdat$sampsize)-length(exdat$sampsize))
fitout<-fitEmaxB(exdat$rslt,exdat$dose,prior,modType=4,prot=exdat$protid,
				count=exdat$sampsize,msSat=msSat,mcmc=mcmc)
plot(fitout)
}
}

\keyword{ dose response}
\keyword{ mcmc }
