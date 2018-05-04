\name{sigmaEmax}
\alias{sigma.fitEmax}
\alias{sigma.fitEmaxB}
\alias{sigma.emaxsim}
\alias{sigma.emaxsimB}

\title{Extract Emax model residual SD estimates}
\description{
Extract Emax model residual SD estimates.  
}
\usage{
	\method{sigma}{fitEmax}(object, ...)
	\method{sigma}{fitEmaxB}(object, ...)
	\method{sigma}{emaxsim}(object, ...)
	\method{sigma}{emaxsimB}(object, ...)
}
\arguments{
	\item{object}{ Output of Emax fitting and simulation functions }
	\item{\dots}{ None additional inputs supported }
}

\value{
	MLE estimate of the residual SD from \code{fitEmax}.  
	Vector of MLE estimates of the residual SD for each \code{emaxsim} simulation.	
	Vector of MCMC generated residual SD for \code{fitEmaxB}.
	Vector of posterior median estimates of the residual SD for 
	each \code{emaxsimB} simulation. 
}

\author{Neal Thomas }

\seealso{ \code{\link{fitEmax}}, \code{\link{fitEmaxB}},
	\code{\link{emaxsim}},\code{\link{emaxsimB}} }

\examples{
doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)
sigma(testout)
}
\keyword{nonlinear}
