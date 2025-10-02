\name{targetCI}
\alias{targetCI}

\title{Compute the dose with confidence interval exceeding a target
change from placebo for each simulated example in an emaxsim object. }
\description{
Selects the lowest dose from a user-specified grid of doses with confidence interval exceeding a targetted
change from placebo for each simulated data set in an emaxsim object.  
}
\usage{
targetCI (object,
	target,
	dgrid,
	clev=0.90,
	high= TRUE)
}
\arguments{
  \item{object}{ An emaxsim object }
  \item{target}{Target improvement from placebo}
  \item{dgrid}{The lowest dose is found by a search over a
user-specified grid of doses.  If dgrid is a single value, it is
interpreted as the number of equally-spaced doses to select from 
zero to the highest dose in the simulated design.}
  \item{clev}{One-sided confidence interval level.}
  \item{high}{When TRUE, lower bounds are computed and must be higher
than the target.  When FALSE, upper bounds must be less than the target.} 
}

\value{

Returns a vector with the lowest dose meeting the criteria.  If a
simulated example does not have a qualifying dose, Inf is returned.
}
\note{
If the grid is very large (>200), execution will slow as a large 
number of estimates and SEs are computed.
}

\author{Neal Thomas }

\seealso{\code{\link{emaxsim}, \link{predict.emaxsim}}, \link{targetD}}
\examples{
	\dontrun{

		# emaxsim changes the random number seed
		nsim<-100
		doselev<-c(0,5,25,50,100)
		n<-c(78,81,81,81,77)

		### population parameters for simulation
		e0<-2.465375 
		ed50<-67.481113 
		
		dtarget<-100
		diftarget<-9.032497
		emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

		sdy<-7.967897
		pop<-c(log(ed50),emax,e0)    
		meanlev<-emaxfun(doselev,pop)  

		###FixedMean is specialized constructor function for emaxsim
		gen.parm<-FixedMean(n,doselev,meanlev,sdy)  

		D1 <- emaxsim(nsim,gen.parm,modType=3)

		target<-6
		tD<- ( (target*ed50)/(emax-target) )
		selectedDose<-targetCI(D1,target,dgrid=c(1:100)+0.5,clev=0.80,high=TRUE)
	}
	\dontshow{

		# emaxsim changes the random number seed
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
		pop<-c(log(ed50),emax,e0)    
		meanlev<-emaxfun(doselev,pop)  

		###FixedMean is specialized constructor function for emaxsim
		gen.parm<-FixedMean(n,doselev,meanlev,sdy)  

		D1 <- emaxsim(nsim,gen.parm,modType=3,nproc=1)

		target<-6
		tD<- ( (target*ed50)/(emax-target) )
		selectedDose<-targetCI(D1,target,dgrid=c(1:100)+0.5,clev=0.80,high=TRUE)
	}
}
\keyword{nonlinear}

