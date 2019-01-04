\name{predict.emaxsim}
\alias{predict.emaxsim}
\title{Mean response and SE for specified doses for each replicate
data set in an emaxsim object }
\description{
Estimated mean/proportion and standard error for each simulated data set in an
emaxsim object.  Also returns mean difference with placebo and their
standard errors.
}
\usage{
\method{predict}{emaxsim}(object, 
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
A list containing:

\item{fitpredv}{Matrix with mean dose response estimate for each
simulated data set.  Number of columns is the number of doses
specified.}
\item{fitdifv}{Matrix with mean dose response estimate minus mean
placebo response for each
simulated data set.  Number of columns is the number of doses
specified.}
\item{sepredv}{Matrix of SEs for \code{fitpredv}.}
\item{sedifv}{Matrix of SEs for \code{fitdifv}.}
}
\author{ Neal Thomas }
\seealso{ \code{\link{emaxsim}}, \code{\link{summary.emaxsim}},
\code{\link{plot.emaxsim}} }
\examples{

\dontrun{
## random number seed changed by this example
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

predout<-predict(D1,c(75,150))
}
\dontshow{
## random number seed changed by this example
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

predout<-predict(D1,c(75,150))
}
}
\keyword{nonlinear}
