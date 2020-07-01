\name{emaxsim}
\alias{emaxsim}

\title{Simulate Emax maximum likelihood estimation}
\description{
 Simulate dose response data and apply 4- or 3- parameter Emax MLE
estimation.  For binary data, the model is fit on the logit scale
and then back-transformed.  When
MLE estimation fails, models with fewer parameters (including models
linear in their parameters) are substituted.  Summaries of estimation
performance are returned for further analyses.  An MCP-MOD test is
also performed for each simulated data set.  
}

\usage{
emaxsim(
	nsim, 
	genObj, 
	modType=3,
	binary=FALSE,
	seed=12357,
	nproc = parallel::detectCores(),
	negEmax=FALSE,
	ed50contr=NULL,
	lambdacontr=NULL,
	testMods=NULL,
	idmax=length(doselev), 
	iparm=NA,
	ed50cutoff=2.5*max(doselev), 
	ed50lowcutoff=doselev[2]/1000,
	switchMod= TRUE,
	truncLambda=6,
	description="")
}

\arguments{
	\item{nsim}{Number of simulation replications}
	\item{genObj}{ Object containing inputs and function to create simulated  
data sets.  These objects are created by special constructor
functions; the current choices are \code{\link{FixedMean}} and \code{\link{RandEmax}}.
}
	\item{modType}{When modType=4, the fitting begins with the 4
parameter model.  If estimation fails or modType=3, the 3-parameter estimation
is applied.  If it fails, a best-fitting model linear in its parameters is
selected.}
	\item{binary}{When specified, the Emax model is fit on the logit scale, and then
 the results are back-transformed to proportions.}
	\item{seed}{Seed for random number generator used to create data.} 
	\item{nproc}{The number of processors to use in parallel computation of the
	simulations, which are divided into equal-sized computational blocks.  When \code{nproc=1}
	a single local processor.}
	\item{negEmax}{When \code{TRUE}, the intended effect is assumed to be negative.}
	\item{ed50contr}{A vector of ED50 values for creating a global null test using the
			  MCP-MOD package DoseFinding based on Emax model-based contrasts.
			  The default is 3 contrasts: the mid-point between pbo and the lowest
			  dose, the mid-point between the 2 highest doses, and the 
			  median of the dose levels.  When there are <=4 doses including pbo,
			  the median-based contrast is excluded.}
	\item{lambdacontr}{Hill parameters matched to the ed50contr.  The default value is 1
for each contrast model.}
	\item{testMods}{The model object for a MCP-MOD test
	created by \code{\link{Mods}} from package \code{DoseFinding}.  If 
specified, the other contrast inputs are ignored.  The  \code{\link{Mods}} 
call should
use the unique sorted dose levels.  The direction of the trend should be specificed in the call to \code{\link{Mods}}.  The \code{negEmax} is stored for use by support functions, but it does not determine the direction of the effect when \code{testMods} is specified.  The validity of
\code{testMods} is not checked.}
	\item{idmax}{Index of the default dose group for comparison to placebo.  
Most analysis functions allow other dose groups to be
specified.  The default is the index of the highest dose.}
	\item{iparm}{Starting values for the  Emax fit.  If
unspecified, starting values are computed.  The
order of the variables is (log(ED50),Emax,E0) or
(log(ED50),lambda,Emax,E0).  Note the transformation of ED50.}
	\item{ed50cutoff}{The upper limit for the ED50 parameter estimates.The default is large
enough to ensure a near linear fit to the data from an 
Emax model.}
	\item{ed50lowcutoff}{Lower allowed limit for the ED50 estimates.}
	\item{switchMod}{If switchMod is TRUE, the algorithm substitutes a
simpler model if (1) convergence is not achieved, (2) the information
matrix is not positive definite at the converged values, (3) the ED50
estimates are outside the cutoff bounds.  If switchMod is F, only
conditions (1) or (2) cause a simpler model to be used.}
  \item{truncLambda}{When \code{modType=4} and the converged estimate of
	  the Hill parameter lambda exceeds \code{truncLambda}, the model fit
	  is judged unstable and discarded.  Set \code{truncLambda=Inf} for 
	  no truncation.  Four parameter model fits are also discarded
	  when lambda is less than 0.1.}
	\item{description}{Optional text describing the simulation setting that
is stored with the simulation output.}
}
\details{
Continuous data can be simulated from any dose response curve with homogeneous
normally distributed residuals. The estimation procedure starts with
ML estimation of a 4- or 3- parameter Emax model depending on modType.
 If modType=3 or 4-parameter estimation fails, a 3 parameter Emax model is fit by
maximum likelihood non-linear least squares.  If 1) nls fails to converge for a 3 parameter Emax model, 2) the
ED50 estimate is <=0, or 3) the ED50 estimate exceeds ed50cutoff, a
linear, log-linear (offset of 1.0), or scaled exponental
(exp(dose/max(dose))), is fit using simple linear least squares
estimation.  The model selected has the smallest residual SS.    

Binary data are handled similarly using maximum likelihood implemented
with the nlm function.  The models are fit on the logit scale and then 
back-transformed for estimation of dose response.  Reduced linear models are 
selected based on the corresponding likelihood deviance.  

MCP-MOD tests are created from contrasts based on the Emax function using 
the \code{DoseFinding} package.  Different
ED50 and lambda (Hill) parameters can be specified to form the contrasts. A contrast
matrix output from the DoseFinding package can be specified instead, allowing for
other contrast choices.  
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
	\item{init}{Matrix with rows containing the starting Emax parameters
for each simulation}
	\item{est4}{Matrix with 4 parmameter Emax fit.  NA if failed to
converge or modType=3}
	\item{est3}{Matrix with 3 parmameter Emax fit.  NA if failed to
converge or 4-parameter estimation was successful.}
	\item{estA}{Matrix with alternative parameter estimates.  NA if
Emax model fit successfully}
	\item{vc}{Variance-covariance matrix for the
estimated parameters stored as a
vector for each simulation.  The vc vector stored has 16,9, or 4
elements 
depending on fitType (with NA values on the end if elements are unused).}
	\item{residSD}{The residual SD based on the selected model.}
	\item{fitType}{Character vector with "4", "3", "L", "LL", or "E" for
4-Emax, 3-Emax, linear,
log-linear, or exponential when an alternative model is selected.}
	\item{pVal}{The \code{nsim} p-values from the global null test.  The p-values are 
	1-sided computed using MCP-Mod.}
	\item{selContrast}{The index of the test contrast producing the smallest p-value.}
	\item{testMods}{Object of class Mods from R package \code{DoseFinding} that defines
	the contrasts used in MCP-MOD testing.  The functions can be plotted with
	\code{DoseFinding} loaded.}
	\item{negEmax}{User input stored for subsequent reference.}
	\item{ed50cutoff}{Upper allowed limit for ED50 estimates}
	\item{ed50lowcutoff}{Lower allowed limit for the ED50 estimates.}
	\item{switchMod}{If switchMod is TRUE, the algorithm substitutes a
simpler model if (1) convergence is not achieved, (2) the information
matrix is not positive definite at the converged values, (3) the ED50
estimates are outside the cutoff bounds.  If switchMod is F, only
conditions (1) or (2) cause a simpler model to be used.}
	\item{negC}{negC=TRUE if the primary fit (from modType) yielded a
ED50 estimate < ED50 lower limit.}
	\item{bigC}{bigC=TRUE if the primary fit (from modType) yielded an 
		  ED50> ED50 upper limit.}
	\item{predpop}{Matrix with population means for each dose group}
	\item{mv}{Matrix with rows containing dose group sample means}
	\item{sdv}{Matrix with rows containing dose group sample SD}
	\item{fitpredv}{Matrix with rows containing dose groups means
estimated from the model}
	\item{sepredv}{Matrix with rows containing SE for fitpredv}
	\item{sedifv}{Matrix with rows containing SE for model-based
differences with placebo}
	\item{rseed}{Starting random
number seed for each simulated data set set that can be assigned to \code{.Random.seed}.
To reproduce the data, the random number generator must also be changed to
\code{RNGkind("L'Ecuyer-CMRG")}.}
	\item{idmax}{Index of default dose group for comparison to placebo
(e.g., for plotting Z-statistics).}
}


\author{ Neal Thomas }

\seealso{  \code{\link{print.emaxsim}}, 
\code{\link{summary.emaxsim}}, \code{\link{plot.emaxsim}}, 
\code{\link{coef.emaxsim}}, \code{\link{sigma.emaxsim}},
\code{\link{vcov.emaxsim}}, \code{\link{predict.emaxsim}},
\code{\link{emaxfun}}}


\examples{

\dontrun{
## emaxsim changes the random number seed
nsim<-50
idmax<-5
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
Ndose<-length(doselev)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 

dtarget<-100
diftarget<-2.464592
emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen<-FixedMean(n,doselev,meanlev,sdy)  

D1 <- emaxsim(nsim,gen,modType=3)
summary(D1,testalph=0.05)

D4 <- emaxsim(nsim,gen,modType=4)
summary(D4,testalph=0.05)
}
\dontshow{
## emaxsim changes the random number seed
nsim<-2
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
Ndose<-length(doselev)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 

dtarget<-100
diftarget<-2.464592
emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen<-FixedMean(n,doselev,meanlev,sdy)  

D1 <- emaxsim(nsim,gen,modType=3,nproc=1)
}
}
\keyword{nonlinear}

