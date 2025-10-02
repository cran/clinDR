\name{emaxsolve}
\alias{solveEmax}
\alias{solveDose}

\title{Solve Emax function for target value}
\description{
Solve the Emax function for dose or Emax to yield a specified response.
}
\usage{
solveEmax(target,dose,led50,lambda,e0,pboadj=TRUE)
solveDose(target,led50,lambda,emax,e0,pboadj=TRUE)
}
\arguments{
	\item{target}{ The targetted response.  If the Emax model is specified
 on the logit scale for binary data, \code{target} and 
 \code{e0} must be logit transformed also.}
	\item{dose}{The dose yielding target.  It is specifed for \code{solveEmax}, and 
			returned for \code{solveDose}}
	\item{led50, lambda, e0}{Emax model parameters (ed50 log transformed)}
	\item{emax}{The Emax model parameter for solveDose.  The value returned
			for \code{solveEmax}}
	\item{pboadj}{When \code{TRUE}, target is placebo-adjusted.}
}

\author{Neal Thomas}

\seealso{\code{\link{fitEmax}}, \code{\link{fitEmaxB}},
	\code{\link{emaxsim}}, \code{\link{emaxsimB}} }

\examples{

e0<-10
dose<-1
led50<-log(0.5)
lambda<-2
target<- -1.5
emax<-solveEmax(target,dose,led50,lambda,e0)
emax

dose1<-solveDose(target,led50,lambda,emax,e0)
dose1

emaxfun(dose=dose1,parm=c(led50,lambda,emax,e0)) - e0

}
\keyword{nonlinear}

