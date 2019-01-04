\name{emaxalt}
\alias{emaxalt}

\title{ Fit 4- or 3-parameter Emax model substituting simpler curves
if convergence not achieved.}
\description{
  ML estimation for 4- and 3-parmeter Emax model.  If the 4-parameter
model is requested, it is estimated and the 3-parameter model is fit
only if  the
4-parameter estimation fails.  If 3-parameter estimation fails, the
linear, log-linear, or exponential model producing the smallest
residual SS is substituted. For binary data, the model 
is fit on the logit scale and then back-transformed.
}
\usage{
emaxalt(y, dose, modType=3,binary=FALSE, 
	iparm=NA,ed50cutoff=2.5*max(doselev), 
	ed50lowcutoff=doselev[2]/1000,switchMod= TRUE,
	truncLambda=6)
}

\arguments{
  \item{y}{Response vector}
  \item{dose}{Doses corresponding to \code{y}}
 \item{modType}{When modType=4, the fitting begins with the 4
parameter model.  If estimation fails or modType=3, the 3-parameter
estimation is applied.  If it fails, a best-fitting model linear in
its parameters is selected.}
	\item{binary}{When specified, the Emax model is fit on the logit scale, and then
 the results are back-transformed to proportions.}
  \item{iparm}{ Vector of optional initial values for the Emax fit. Starting 
  values are computed if not specified. }
  \item{ed50cutoff}{ Upper allowed limit for ED50 estimates. }
  \item{ed50lowcutoff}{Lower allowed limit for the ED50 estimates.}
  \item{switchMod}{If switchMod is TRUE, the algorithm substitutes a
simpler model if (1) convergence is not achieved, (2) the information
matrix is not positive definite at the converged values, (3) the ED50
estimates are outside the cutoff bounds.  If switchMod is F, only
conditions (1) or (2) cause a simpler model to be used.}
  \item{truncLambda}{When \code{modType=4} and the converged estimate of
	  the Hill parameter lambda exceeds \code{truncLambda}, the model fit
	  is judged unstable and discarded.  Set \code{truncLambda=Inf} for 
	  no truncation.}
}
\details{
  The partial linear method is used in nls.  If it fails, gauss-newton
is attempted.  If both methods fail, the next simpler model is
attempted.  For the 4-parameter model, the next step is the
3-parameter model. For the 3-parameter model, a linear, log-linear
log(dose+1.0), and exp(dose/max(dose)) are fit using lm, and the
2-parm fit with the smallest residual SS is selected.  
}
\value{
  A list assigned class "emaxalt" with the following elements:

  \item{dm}{Vector containing dose group means}
  \item{dsd}{Vector containing dose group SDs}
  \item{Sparm}{Vector of starting values for 3-parameter Emax fit.}
  \item{fitType}{Character vector with "4", "3", "L", "LL", or "E" for
4-Emax, 3-Emax, linear,
log-linear, or exponential when an alternative model is selected.}
  \item{vc}{The variance-covariance matrix of the model parameters
stored as a vector.  The length is 16, 9, 4 depending on fitType.}
  \item{fitpred}{Dose groups means
estimated from the model}
  \item{residSD}{The residual SD based on the selected model.}
 \item{sepred}{SEs for estimates in
fitpred}
  \item{sedif}{SEs for model-based estimates of difference with
placebo}
  \item{bigC}{bigC= TRUE if the primary fit (from modType) yielded an 
              ED50 >ED50 upper limit.}
  \item{negC}{negC= TRUE if the primary fit (from modType) yielded a
ED50 estimate < ED50 lower limit.}
  \item{est4}{ 4-parmameter Emax fit (ed50,lambda,emax,e0).  NA if
failed to converge or 3-parameter model requested.}
  \item{est3}{ 3-parmameter Emax fit (ed50,emax,e0).  NA if
failed to converge or 4-parameter model successfully fit.}
 \item{estA}{Alternative parameter estimates.  NA if
Emax model fit successfully}
}

\author{Neal Thomas }

\seealso{ \code{\link{emaxsim}}, \code{\link{nls} }}
\examples{

save.seed<-.Random.seed
set.seed(12357)

doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
dose<-rep(doselev,n)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanresp<-emaxfun(dose,pop)  
y<-rnorm(sum(n),meanresp,sdy)

simout<-emaxalt(y,dose)

simout2<-emaxalt(y,dose,modType=4)

.Random.seed<-save.seed
}
\keyword{nonlinear}

