\name{print.emaxsim}
\alias{print.emaxsim}
\title{Print simulation output from emaxsim }
\description{
  Prints key summary variables of Emax estimation peformance for each
simulation.  Can be used to identify simulated data sets yielding
problems with common estimation methods.
}
\usage{
\method{print}{emaxsim}(x, 
        nprint = min(length(x$fitType), 20), 
        id = x$idmax, 
        digits = 3, ...)
}

\arguments{
  \item{x}{Output of \code{\link{emaxsim}} }
  \item{nprint}{ Number of simulations to print.  If a vector of
length 2, \code{nprint} is the range of simulations to print. }
  \item{id}{Output includes the stdBias for the dose with index \code{id}  vs
placebo }
  \item{digits}{ Number of decimal digits to print for Z and p-values }
  \item{\dots}{ Other print parameters (none currently implemented) }
}

\value{
Printed output returned as invisible matrix.

}
\author{ Neal Thomas }
\note{ 
The stdBias printed is the difference between the estimated
dose response at the dose with index \code{id} and its population
value.  The difference is divided by the SE of the estimator computed
using the delta method.
}
\seealso{ \code{\link{emaxsim}}, \code{\link{summary.emaxsim}},
\code{\link{plot.emaxsim}} }
\examples{


\dontrun{
## emaxsim changes the random number seed
nsim<-50
idmax<-5
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop.parm<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop.parm)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen.parm)

print(D1,c(31,50),digits=2,id=4)

print(D1,c(1,20))

D1  ### implicitly calls print with default parameter settings
}
\dontshow{
## emaxsim changes the random number seed
nsim<-3
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop.parm<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop.parm)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen.parm,nproc=1)

print(D1,c(2,3),digits=2,id=4)
}
}
\keyword{nonlinear}
