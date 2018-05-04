\name{vcovEmax}
\alias{vcov.fitEmax}
\alias{vcov.emaxsim}

\title{Extract Emax model variance-covariance matrix for ML estimates}
\description{
Extract Emax model variance-covariance matrix for ML estimates
}
\usage{
	\method{vcov}{fitEmax}(object, ...)
	\method{vcov}{emaxsim}(object, ...)
}
\arguments{
	\item{object}{ Output of Emax fitting and simulation functions }
	\item{\dots}{ None additional inputs supported }
}

\value{
	Variance-Covariance matrix for the MLE estimates of the parameters from \code{fitEmax}. 
	The lower half of the variance-covariance matrix for the estimated 
	parameters stored as a vector in column-major order for each 
	\code{emaxsim} simulation. The vc matrix has 16,9, or 4 
	elements depending on fitType. 
}

\author{Neal Thomas }

\seealso{ \code{\link{fitEmax}}, 
	\code{\link{emaxsim}} }

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
vcov(testout)
}
\keyword{nonlinear}
