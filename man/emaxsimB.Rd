\name{emaxsimB}
\alias{emaxsimB}
\title{
Simulate Emax Bayesian estimation
}
\description{

Simulate dose response data and apply 4- or 3- parameter
sigmoidal or hyperbolic Bayesian estimation.  The prior distributions
are input by the user based on the empirical distributions estimated
from dose response meta-analyses.  For binary response data, the Emax
model is fit on the logit scale, and then back-transformed}
\usage{
emaxsimB(nsim, genObj, prior, modType = 3, binary = FALSE, seed=12357,
	check = FALSE, nproc=parallel::detectCores(), negEmax = FALSE, ed50contr = NULL, 
	lambdacontr = NULL, testMods = NULL, idmax = length(doselev), 
	mcmc = mcmc.control(), 
	customCode=NULL, customParms=NULL,
	description = "")
}
\arguments{
  \item{nsim}{
Number of simulation replications
}
  \item{genObj}{
Object containing inputs and function to create simulated  
data sets.  These objects are created by special constructor
functions; the current choices are \code{\link{FixedMean}} and \code{\link{RandEmax}}.
}
  \item{prior}{
Prior distribution input using \code{\link{prior.control}}
}
  \item{modType}{
Specifies sigmoidal or hyperbolic emax model.
}
  \item{binary}{
When specified, the Emax model is fit on the logit scale, and then
 the results are back-transformed to proportions.
}
	\item{seed}{Seed for random number generator used to create data.  A separate
	seed can be passed to \code{rstan} through the \code{MCMC} object.}
  \item{check}{
When \code{TRUE}, a single simulated data set is created and the data
and \code{rstan} object are returned for convergence checking.  The data
are in the form needed for developing \code{customCode}.  Note that
\code{customCode} is not called when \code{check=TRUE}.
}
	\item{nproc}{The number of processors to use in parallel computation of the
	simulations, which are divided into equal-sized computational blocks.  When \code{nproc=1}
	a single local processor.}
  \item{negEmax}{
When \code{TRUE}, the intended effect is assumed to be negative.
}
	\item{ed50contr}{A vector of ED50 values for creating a global null test using the
			  MCP-MOD package DoseFinding based on Emax model-based contrasts.
			  The default is 3 contrasts: the mid-point between pbo and the lowest
			  dose, the mid-point between the 2 highest doses, and the 
			  median of the dose levels.  When there are <=4 doses including pbo,
			  the median-based contrast is excluded.}
	\item{lambdacontr}{Hill parameters matched to the ed50contr.  The default value is 1
for each contrast model.}
	\item{testMods}{The model object for a MCP-MOD test
	created by \code{\link{Mods}} from package DoseFinding.  If 
specified, the other contrast inputs are ignored.  The  \code{\link{Mods}} 
call should
use the unique sorted dose levels.  The direction of the trend should be specificed in the call to \code{\link{Mods}}.  The \code{negEmax} is stored for use by support functions, but it does not determine the direction of the effect when \code{testMods} is specified.  The validity of
\code{testMods} is not checked.}
	\item{idmax}{Index of the default dose group for comparison to placebo.  
Most analysis functions allow other dose groups to be
specified.  The default is the index of the highest dose.}
  \item{mcmc}{
MCMC settings created using \code{\link{mcmc.control}}
}
	\item{customCode}{An optional user supplied function that computes custom
	estimates/decision criteria from each simulated data set and its Bayesian
	model fit.  The output are stored in a list, \code{customOut}, of length
	\code{nsim}.  See the Details section below for a description of the
	mandatory inputs to the \code{customCode} function.
	}
	\item{customParms}{Optional parameters that can be passed to 
	\code{customCode}.}
  \item{description}{
  Optional text describing the simulation setting that
is stored with the simulation output.
}

}
\details{
The Bayesian model fits are implemented in \code{rstan} using function
\code{\link{fitEmaxB}}.  The function \code{compileStanModels} must be 
executed once to create compiled STAN
code before \code{emaxsimB} can be used.

Continuous data can be simulated from any dose response curve with homogeneous
normally distributed residuals.

Binary data are handled similarly.  The models are fit on the logit scale and then 
back-transformed for estimation of dose response.  Reduced linear models are 
selected based on the corresponding likelihood deviance.  

MCP-MOD tests are created from contrasts based on the Emax function using
the \code{DoseFinding} package.  Different
ED50 and lambda (Hill) parameters can be specified to form the contrasts. A contrast
matrix output from the DoseFinding package can be specified instead, allowing for
other contrast choices.  

Customized code:

For binary data, the inputs to the function customCode for each simulated data set
will be (parms,pVal,dose,y), where parms is the matrix of parameters
generated from the posterior distribution with columns in the order given in
function \code{\link{emaxfun}}, pVal is the MCP-MOD p-value, dose and y are
the patient-level simulated data.  For continuous data, the inputs
are (parms,residSD,pVal,dose,y), where \code{residSD} are the variance
parameters generated from their posterior distribution.  The \code{customParms} supply other user-inputs 
such as a target efficacy level.  When it is not null, the \code{customCode}
inputs must be (parms,pVal,dose,y,customParms) or (parms,residSD,pVal,dose,y,customParms).
}
\value{
  A list is returned with class(emaxsim) containing:

	\item{description}{User description of simulation}
	\item{binary}{Binary response data.}
	\item{modType}{User supplied starting Emax model}
	\item{genObj}{List object with data and function used to generate study
data}
	\item{pop}{Matrix with rows containing population parameters
for each simulation.  Type of parameter depends on
constructor function generating study data.}
	\item{popSD}{Vector containing the population SD used to generate
	continuous data.  \code{NULL} for binary data.}
	\item{mcmc}{mcmc input settings}
	\item{prior}{Input prior distribution.}
	\item{est}{Matrix with posterior median parameter estimates for each simulation:
	(led50,lambda,emax,e0) or (led50,emax,e0).}
	\item{residSD}{The posterior median of the residual SD for each simulation.}
	\item{pVal}{The \code{nsim} p-values from the global null test.  The p-values are 
	1-sided computed using MCP-Mod.}
	\item{selContrast}{The index of the test contrast producing the smallest p-value.}
	\item{testMods}{Object of class Mods from R package \code{DoseFinding} that defines
	the contrasts used in MCP-MOD testing.  The functions can be plotted with
	\code{DoseFinding} loaded.}	
	\item{gofP}{Goodness of fit test computed by \code{checkMonoEmax}.}
	\item{negEmax}{User input stored for subsequent reference.}
	\item{predpop}{Matrix with population means for each dose group}
	\item{mv}{Matrix with rows containing dose group sample means}
	\item{sdv}{Matrix with rows containing dose group sample SD}
	\item{msSat}{Pooled within-dose group sample variance}
	\item{fitpredv}{Matrix with rows containing dose groups means
estimated by the posterior means.}
	\item{sepredv}{Matrix with rows containing SE (posterior SD) associated with fitpredv}
	\item{sedifv}{Matrix with rows containing SE (posterior SD) for the posterior means of the
differences with placebo}
	\item{lb,ub}{Array with lower posterior (0.025,0.05,0.1) and upper posterior
(0.975,0.95,0.9) percentiles of differences between dose group means and 
placebo.  The array ordering is dose group minus placebo, simulation, and percentile.}
	\item{rseed}{Starting random
number seed for each simulated data set set that can be assigned to \code{.Random.seed}.
To reproduce the data, the random number generator must also be changed to
\code{RNGkind("L'Ecuyer-CMRG")}.}
	\item{idmax}{Index of default dose group for comparison to placebo
(e.g., for plotting Z-statistics).}
	\item{customOut}{List with customized output. It will be \code{NULL} if \code{customCOde} is not specified.}
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


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{  \code{\link{print.emaxsimB}}, 
\code{\link{summary.emaxsimB}}, \code{\link{plot.emaxsimB}}, 
\code{\link{coef.emaxsimB}}, \code{\link{sigma.emaxsimB}},
\code{\link{emaxfun}}}

\examples{
\dontrun{

### emaxsimB changes the random number seed

nsim<-50
idmax<-5
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
Ndose<-length(doselev)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-4.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen<-FixedMean(n,doselev,meanlev,sdy)  

prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
		sigmaup=30,edDF=5)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,
		propInit=0.15,adapt_delta = 0.95)

### custom code to compute the distribution of the dose yielding
### a target diff with pbo
customCode<-function(parms,residSD,pVal,dose,y,customParms){
	target<-customParms
	ed50<-exp(parms[,1])
	emax<-parms[,2]
	td<-ifelse(emax-target>0,ed50/(emax-target),Inf)
	tdest<-median(td)
	lb<-quantile(td,0.1)
	ub<-quantile(td,0.9)
	return(c(td=tdest,lb=lb,ub=ub))
}

D1 <- emaxsimB(nsim,gen, prior, modType=3,seed=12357,mcmc=mcmc,check=FALSE,
				customCode=customCode,customParms=1.0)
D1
}
}
\keyword{ Bayes}
\keyword{ Emax }
