\name{plot.emaxsim}
\alias{plot.emaxsim}

\title{  Plot the output of emaxsim}
\description{
  A Q-Q plot of the dose response estimate of the mean at a specified dose minus
the population value divided by the standard error of the estimator
(computed using the delta method)
}
\usage{
\method{plot}{emaxsim}(x, id = x$idmax, difPlacebo= TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Output  of \code{\link{emaxsim}} }
  \item{id}{ Index of the dose to be assessed (placebo index=1). }
  \item{difPlacebo}{ If true (default), the estimates and population
values are differences with placebo.  IF false, absolute dose response
values are used.}
  \item{\dots}{ Optional parameters passed to the plotting function }
}

\value{
No output is returned.
}

\author{Neal Thomas }

\seealso{  \code{\link{emaxsim}}, \code{\link{print.emaxsim}},
\code{\link{summary.emaxsim}}  }

\examples{

\dontrun{
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

plot(D1,id=3)
}
\dontshow{
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
D1 <- emaxsim(nsim=3,gen.parm,nproc=1)

plot(D1,id=3)
}
}
\keyword{nonlinear}
