\name{predict.emaxsimobj}
\alias{predict.emaxsimobj}
\title{Mean response and SE for specified doses for a
simulated emaxsimobj object}
\description{
Estimated mean/proportion and standard error for specified doses computed from
the output of a simulated data set created by function emaxsim.  Also returns mean difference with placebo and their
standard errors.
}
\usage{
\method{predict}{emaxsimobj}(object, 
          dose, dref=0,
          ...)
}

\arguments{
  \item{object}{Output of the extract function [] applied to an object
  createad by \code{\link{emaxsim}}. }
  \item{dose}{ Vector (can be a single value) of doses where dose
response curve is to be evaluated. }
   \item{dref}{A reference dose (0 by default) for contrasts, but other values 
   can be specified.  If specified, a single reference value must be given.
   }
  \item{...}{Optional arguments are not used.}
}

\value{
A list containing:

\item{fitpred}{Vector with mean dose response estimate for each
specified dose.}
\item{fitdif}{Corresponding differences with placebo.}
\item{sepred}{ SEs for \code{fitpred}.}
\item{sedif}{SEs for \code{fitdif}.}
}
\author{ Neal Thomas }
\seealso{ \code{\link{emaxsim}}, \code{\link{summary.emaxsim}},
\code{\link{predict.emaxsim}} }
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

dtarget<-100
diftarget<-9.032497
emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

sdy<-7.967897
pop.parm<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop.parm)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen.parm)
d10<-D1[10]
predict(d10,c(75,150))
}
\dontshow{
## emaxsim changes the random number seed
nsim<-3
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 

dtarget<-100
diftarget<-9.032497
emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

sdy<-7.967897
pop.parm<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop.parm)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen.parm,nproc=1)
d3<-D1[3]
predict(d3,c(75,150))
}
}
\keyword{nonlinear}



