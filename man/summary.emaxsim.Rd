\name{summary.emaxsim}
\alias{summary.emaxsim}

\title{ Summary of output of emaxsim }
\description{
Detailed summary of repeated sampling properties of Emax estimation
and comparison with simple pairwise comparisons.
}
\usage{
\method{summary}{emaxsim}(object, testalpha = 0.05, clev = 0.9, 
                          seSim = FALSE, ...)
}

\arguments{
  \item{object}{ Output  of \code{\link{emaxsim}} }
  \item{testalpha}{Alpha level for a one-sided MCP-MOD trend test}
  \item{clev}{ Nominal confidence level for reported CIs }
  \item{seSim}{ If \code{TRUE}, then simulation standard errors are reported in
parentheses.  These should be distinguished from standard errors for
estimators in the simulation.}
  \item{\dots}{Other unspecified parameters (none currently utilized) }
}

\value{
The function produces annotated output summarizing the
repeated sampling properties of the estimation procedures.  It does not return any output.
}
\author{ Neal Thomas }

\seealso{  \code{\link{emaxsim}}, \code{\link{print.emaxsim}},
\code{\link{plot.emaxsim}}  }

\examples{

\dontrun{
## emaxsim changes the random number seed
nsim<-50
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
summary(D1,testalph=0.05,clev=0.95,seSim=TRUE)
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
summary(D1,testalph=0.05,clev=0.95,seSim=TRUE)
}
}
\keyword{nonlinear}
