\name{predict.emaxalt}
\alias{predict.emaxalt}
\title{Mean response and SE for specified doses for a
simulated object output by function emaxalt}
\description{
Estimated mean and standard error for specified doses computed from
the output of a model fit by function
emaxalt.  Also returns mean difference with placebo and their
standard errors.
}
\usage{
\method{predict}{emaxalt}(object,dose, dref=0, ...)
}

\arguments{
  \item{object}{Output of \code{\link{emaxalt}} }
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
\seealso{ \code{\link{emaxalt}}, \code{\link{predict.emaxsimobj}},
\code{\link{predict.emaxsim}} }
\examples{

\dontrun{
## random number seed changed by this example

doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
dose<-rep(doselev,n)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop.parm<-c(log(ed50),e0,emax)    
meanresp<-emaxfun(dose,pop.parm)  
y<-rnorm(sum(n),meanresp,sdy)

simout<-emaxalt(y,dose)
predict(simout,c(75,150))

simout2<-emaxalt(y,dose,modType=4)
predict(simout2,c(75,150))
}
\dontshow{
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
dose<-rep(doselev,n)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop.parm<-c(log(ed50),e0,emax)    
meanresp<-emaxfun(dose,pop.parm)  
y<-rnorm(sum(n),meanresp,sdy)

simout<-emaxalt(y,dose)
predict(simout,c(75,150))
}
}
\keyword{nonlinear}
