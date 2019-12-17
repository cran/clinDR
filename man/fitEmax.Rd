\name{fitEmax}
\alias{fitEmax}

\title{ML fit of hyperbolic or sigmoidal Emax models to continuous/binary dose
response data.
}
\description{
Calls Newton-Raphson optimizers, nls and nlm,
for a hyperbolic or sigmoidal 
Emax model.  Different intercepts for multiple protocol-data are supported.  For binary
data, the Emax model is on the logit scale.
}
\usage{
fitEmax(y,dose,iparm,xparm,modType=3,
	prot=rep(1,length(y)),count=rep(1,length(y)),xbase=NULL,
	binary=FALSE,diagnostics=TRUE,msSat=NULL,
	pboAdj=rep(FALSE,max(prot)),optObj=TRUE)
}

\arguments{

  \item{y}{Outcome for each patient.  Missing \code{Y} values
are are not permitted.  Dose/protocol group means for grouped continuous data.  For binary
data, \code{y} must be 0/1 and counts must be supplied for each 0/1 value.}
  \item{dose}{Dose for each patient.}
  \item{iparm}{Optional starting values for the Newton-Raphson algorithm.  The
order of the variables is (log(ED50),Emax,E0) or
(log(ED50),lambda,Emax,E0).  Note the transformation of ED50. If there is more
than one protocol, the E0 is automatically duplicated.}
	\item{xparm}{Optional starting values for the baseline covariate slopes (if any).
\code{xparm} must be specified when \code{iparm} and \code{xbase} are specified.
\code{startEmax} is used to obtain starting values if no starting values are specified.}
  \item{modType}{modType=3 (default) for the 3-parameter hyperbolic Emax model.
  modType=4 for the 4-parameter sigmoidal Emax model.} 
  \item{prot}{Protocol (group) membership used to create multiple intercepts.  
  The default is a single protocol.}
  \item{count}{Counts for the number of patients when the \code{Y} are dose 
  continuous group means or binary 0/1 values.  
  Default is 1 (ungrouped data).}
  \item{xbase}{A matrix of baseline covariates with rows corresponding to 
   \code{y} that enter as linear additive predictors.  The baseline covariates must 
  be centered about
  their (protocol-specific) means.  \code{xbase} does not include an intercept or
  protocol indicators.  Covariates cannnot be specified with PBO adjusted or 
  aggregated input.}
  \item{diagnostics}{Print trace information per iteration and any error messages from the optimizing
	methods.  Printing can be suppressed for use in simulation studies.}
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
  \item{optObj}{Include the output object from the R optimization code in the \code{fitEmax}
  output.}
}

\details{
Fits the 3- or 4- Emax model  using
\code{\link{nls}}.  A newton-raphson algorithm is tried first 
followed by a partial linear optimatization if needed.  Binary
data are fit using \code{\link{nlm}}.}

\value{
A list assigned class "fitEmax" with:
\item{fit}{The parameter estimates and their variance-covariance matrix.}
\item{y, dose, modType, prot, count, binary, pboAdj}{Input values.}
\item{gofTest}{Goodness of fit p-value based on likelihood ratio comparison of
	the model to a saturated fit.}
\item{nll}{\code{-2*loglikelihood} for the Emax model and the saturated model.  
	Residual sums of squares are returned for continuous data models.  These
statistics can be used to construct other tests using multiple calls to
\code{fitEmax} (e.g., 3 vs 4 parameter Emax
models, or a common intercept model across protocols).}
\item{df}{Residual degrees of freedom for the Emax model and the saturated model.}
\item{optobj}{When requested, the fit object returned by the R optimation functions.}
}

\author{Neal Thomas }

\seealso{\code{\link{nls}}, \code{\link{nlm}},  \code{\link{nllogis}}, 
\code{\link{predict.fitEmax}}, \code{\link{plot.fitEmax}}, \code{\link{coef.fitEmax}} }
\examples{

## the example changes the random number seed

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)
}

\keyword{nonlinear}
