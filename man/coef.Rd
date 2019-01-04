\name{coefEmax}
\alias{coef.fitEmax}
\alias{coef.fitEmaxB}
\alias{coef.emaxsim}
\alias{coef.emaxsimB}

\title{Extract Emax model parameter estimates}
\description{
Extract Emax model parameter estimates.  MLE for fitEmax. 
Matrix of MCMC generated parameters for fitEmaxB.
}
\usage{
	\method{coef}{fitEmax}(object, ...)
	\method{coef}{fitEmaxB}(object, ...)
	\method{coef}{emaxsim}(object, ...)
	\method{coef}{emaxsimB}(object, ...)
}
\arguments{
	\item{object}{ Output of Emax fitting function }
	\item{\dots}{ None additional inputs supported }
}

\value{
	Vector of MLE estimates of model parameter from \code{fitEmax}.  
	Matrix of MCMC generated parameters for \code{fitEmaxB}.
	Matrix with posterior median parameter estimates for 
	each \code{emaxsimB} simulation: (led50,lambda,emax,e0) or (led50,emax,e0).
	For \code{emaxsim}, a list is returned with the model type fit for each
	simulation, and a matrix with the corresponding model coefficients.  The 
	order of the parameters is given in the \code{emaxsim} documentation.
}

\author{Neal Thomas }

\seealso{ \code{\link{sigma}}, \code{\link{fitEmax}}, \code{\link{fitEmaxB}},
	\code{\link{emaxsim}}, \code{\link{emaxsimB}} }

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
coef(testout)
}
\keyword{nonlinear}
