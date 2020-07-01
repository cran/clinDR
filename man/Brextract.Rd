\name{"Extract.emaxsim"}
\alias{"Extract.emaxsim"}
\alias{[.emaxsim}

\title{Extract a simulation from the output of emaxsim }
\description{
  Extract a simulated data set from the output of
emaxsim. Data are re-created using the stored random number seed.
}

\usage{
\method{[}{emaxsim}(x, i, ...)
}

\arguments{
  \item{x}{Output object from \code{\link{emaxsim}} }
  \item{i}{Simulation replication to extract }
  \item{\dots}{Parameters passed to other functions (none currently) }
}
\details{
Re-creates the ith simulated data set for subsequent analyses.  Also returns all
analyses done for the ith data set in  \code{\link{emaxsim}}
}
\value{
A list is returned with class(emaxsimobj) containing:

  \item{y}{Response vector}
  \item{dose}{Doses corresponding to \code{y}}
  \item{pop}{Population parameters; type of parameter depends on
constructor function generating study data.}
	\item{popSD}{Vector containing the population SD used to generate
	continuous data.  \code{NULL} for binary data.}
  \item{init}{Starting Emax parameters}
  \item{est4}{ 4-parmameter Emax fit (ed50,lambda,emax,e0).  NA if
failed to converge or 3-parameter model requested.}
  \item{est3}{ 3-parmameter Emax fit (ed50,emax,e0).  NA if
failed to converge or 4-parameter model successfully fit.}
 \item{estA}{Alternative parameter estimates.  NA if
Emax model fit successfully}
  \item{vc}{The variance-covariance matrix of the model parameters 
for the selected model. }
  \item{residSD}{The residual SD based on the selected model.}
  \item{bigC}{bigC= TRUE if the primary fit (from modType) yielded an 
              ED50 > ED50 upper limit.}
  \item{negC}{negC= TRUE if the primary fit (from modType) yielded a
negative ED50 estimate< ED50 lower limit }
 \item{modType}{When modType=4, the fitting begins with the 4
parameter model.  If estimation fails or modType=3, the 3-parameter
estimation is applied.  If it fails, a best-fitting model linear in
its parameters is selected.}
  \item{fit}{Output of model determined by fitType}
  \item{fitType}{Character vector with "4", "3", "L", "LL", or "E" for
4-Emax, 3-Emax, linear,
log-linear, or exponential when an alternative model is selected.}
  \item{ed50cutoff}{Upper allowed limit for ED50 estimates.}
  \item{ed50lowcutoff}{Lower allowed limit for the ED50 estimates.}
  \item{switchMod}{If switchMod is TRUE, the algorithm substitutes a
simpler model if (1) convergence is not achieved, (2) the information
matrix is not positive definite at the converged values, (3) the ED50
estimates are outside the cutoff bounds.  If switchMod is F, only
conditions (1) or (2) cause a simpler model to be used.}
   \item{PL}{T if the 'plinear' algorithm in
nls converged}
  \item{predpop}{Population means for each dose group}
  \item{dm}{Vector containing dose group means}
  \item{dsd}{Vector containing dose group SDs}
  \item{fitpred}{Dose groups means
estimated from the model}
 \item{sepred}{SEs for estimates in
fitpred}
  \item{sedif}{SEs for model-based estimates of difference with placebo}
  \item{pVal, selContrast}{P-value and contrast selected from MCP-MOD test}
  \item{idmax}{Index of default dose group for comparison to placebo}
}

\author{Neal Thomas }

\seealso{\code{\link{emaxsim}}, \code{\link{print.emaxsimobj}}, 
\code{\link{plot.emaxsimobj}}, \code{\link{update.emaxsimobj}} }
\examples{

\dontrun{
## code change random number seed

nsim<-50
idmax<-5
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 

dtarget<-100
diftarget<-9.032497
emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  

D1 <- emaxsim(nsim,gen.parm,modType=3)
e49<-D1[49]                  #### extract 49th simulation

}
\dontshow{
## code change random number seed

doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 

dtarget<-100
diftarget<-9.032497
emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  

D1 <- emaxsim(nsim=2,gen.parm,modType=3,nproc=1)
e49<-D1[2]                  
}
}

\note{
Extraction from a simulation object requires re-creation of the
simulated data set.  If the
extracted object is to be used more than once, it is more efficient to
save the extracted object than reuse [].
}
\keyword{nonlinear}

