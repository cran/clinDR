\name{predict.emaxsimBobj}
\alias{predict.emaxsimBobj}
\title{Mean response estimates (posterior means) and SE (posterior SD) 
for specified doses for a simulated emaxsimBobj object}
\description{
Estimated mean and standard error for specified doses (posterior means and SD) computed from
the output of a simulated data set created by function emaxsimB.  
Also returns mean difference with placebo and their
standard errors.
}
\usage{
\method{predict}{emaxsimBobj}(object, 
                dose, dref=0, clev=0.9, 
                ...)
}

\arguments{
  \item{object}{Output of the extract function [] applied to an object
  createad by \code{\link{emaxsimB}}. }
  \item{dose}{ Vector (can be a single value) of doses where dose
response curve is to be evaluated. }
   \item{dref}{A reference dose (0 by default) for contrasts, but other values 
   can be specified.  If specified, a single reference value must be given.
   }
   \item{clev}{Specified probablity of the posterior interval}
  \item{...}{Optional arguments are not used.}
}

\value{
A list containing:

\item{pred}{Vector with mean dose response estimates for each
specified dose.}
\item{fitdif}{Corresponding differences with placebo.}
\item{se}{ SEs (posterior SD) for \code{pred}.}
\item{sedif}{SEs (posterior SD) for \code{fitdif}.}
\item{lb, ub, lbdif, ubdif}{Bounds of \code{clev} posterior intervals.}
}
\author{ Neal Thomas }
\seealso{ \code{\link{emaxsim}}, \code{\link{summary.emaxsim}},
\code{\link{predict.emaxsim}} }
\examples{

\dontrun{
### emaxsimB changes the random number seed
nsim<-50
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
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)

D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)
predict(D1[1],dose=c(75,125))
}
}
\keyword{nonlinear}
