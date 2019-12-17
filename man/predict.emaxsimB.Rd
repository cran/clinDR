\name{predict.emaxsimB}
\alias{predict.emaxsimB}
\title{Mean response and SE for each replicate
data set in an emaxsimB object }
\description{
Return warning and explanation that only predicted values at doses
included in the study are available.   The code needed to obtain
predicted values at other doses is indicated.
}
\usage{
\method{predict}{emaxsimB}(object, 
          dose, dref=0, ...)
}

\arguments{
  \item{object}{Output of \code{\link{emaxsim}} }
  \item{dose}{ Vector (can be a single value) of doses where dose
response curve is to be evaluated. }
   \item{dref}{A reference dose (0 by default) for contrasts, but other values 
   can be specified.  If specified, a single reference value must be given.
   }
  \item{...}{Optional arguments are not used.}
}

\value{
No output.  
}
\author{ Neal Thomas }
\seealso{ \code{\link{emaxsimB}}, \code{\link{summary.emaxsimB}},
\code{\link{plot.emaxsimB}} }
\examples{

\dontrun{
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
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,
		propInit=0.15,adapt_delta = 0.95)


D1 <- emaxsimB(nsim,gen, prior, modType=3,seed=12357,mcmc=mcmc,check=FALSE)

predict(D1,dose=20)
}
}
\keyword{nonlinear}
