\name{"Extract.emaxsimB"}
\alias{"Extract.emaxsimB"}
\alias{[.emaxsimB}

\title{Extract a simulation from the output of emaxsimB }
\description{
  Extract a simulated data set from the output of
emaxsimB. Data are re-created using the stored random number seed.
}

\usage{
\method{[}{emaxsimB}(x, i, ...)
}

\arguments{
  \item{x}{Output object from \code{\link{emaxsimB}} }
  \item{i}{Simulation replication to extract }
  \item{\dots}{Parameters passed to other functions (none currently) }
}
\details{
Re-creates the ith simulated data set for subsequent analyses.  Also returns all
analyses done for the ith data set in  \code{\link{emaxsimB}}
}
\value{
A list is returned with class(emaxsimBobj) containing:

  \item{y}{Response vector}
  \item{dose}{Doses corresponding to \code{y}}
  \item{pop}{Population parameters; type of parameter depends on
constructor function generating study data.}
	\item{popSD}{Vector containing the population SD used to generate
	continuous data.  \code{NULL} for binary data.}
	\item{binary}{When \code{TRUE}, binary data modeled on the logit scale}
 \item{modType}{\code{modType=3, 4}, for the hyperbolic and sigmoidal 
 Emax models.}
	\item{predpop}{Population means for each dose group}
  \item{dm}{Vector containing dose group means}
  \item{dsd}{Vector containing dose group SDs}
  \item{fitpred}{Posterior means of the dose groups means}
	\item{sepred}{SE (posterior SD) corresponding to the estmates in
fitpred}
  \item{sedif}{SE (posterior SD) for the differences with placebo}
  \item{bfit}{Bayesian fitted model of class \code{fitEmaxB}.}
  \item{prior, mcmc}{See \code{fitEmax} for documentation.}
  \item{pVal, selContrast}{P-value and contrast selected from MCP-MOD test}
  \item{idmax}{Index of default dose group for comparison to placebo}
}

\author{Neal Thomas }

\seealso{\code{\link{emaxsimB}}, \code{\link{print.emaxsimBobj}}, 
\code{\link{plot.emaxsimBobj}} }
\examples{

\dontrun{

save.seed<-.Random.seed
set.seed(12357)

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

prior<-emaxPrior.control(epmu=0,epsca=30,difTargetmu=0,
		difTargetsca=30,dTarget=100,p50=50,sigmalow=0.1,
		sigmaup=30,parmDF=5)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)

D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)

out<-D1[2]


.Random.seed<-save.seed
}
}

\note{
Extraction from a simulation object requires re-creation of the
simulated data set.  If the
extracted object is to be used more than once, it is more efficient to
save the extracted object than reuse [].
}
\keyword{nonlinear}

