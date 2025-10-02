\name{targetD}
\alias{targetD}

\title{Compute the MLE (and its SE) of the dose achieving a specified
target improvement from placebo. }
\description{
The MLE (se) of the dose required to achieve a targetted improvement
from placebo.
The fit can be from a 3- or 4- parameter Emax model or output from 
function emaxalt, or an object of class emaxsimobj.  The Emax model
is on the logit scale for binary data.
}
\usage{
targetD (fit,
	target,
	modType=4,
	binary=FALSE)
}
\arguments{
  \item{fit}{Output of \code{\link{nls}} \code{fit} to a 3- or 4-parameter Emax
model.  The order of the parameters in the fit must be (log(ed50),emax,e0)
or (log(ed50),lambda,emax,e0).  
\code{fit} can also be a list with the first element the coefficient
vector, and the second element the variance-covariance matrix. 
Alternatively, \code{fit} may be of class \code{emaxalt} or 
\code{emaxsimobj}, and the target dose is based on the fitted model.}
  \item{target}{Targetted change from placebo (positive or negative). }
  \item{modType}{Value is 3 or 4 for the 3 or 4-parameter Emax model
output from nls with parameters in the
order (ed50,emax,e0) or (ed50,lambda,emax,e0).  
\code{modType} is ignored if fit is from \code{emaxalt} or \code{emaxsimobj}.}
  \item{binary}{When \code{TRUE}, the fit is assumed to be for binary data on 
  the logistic scale.  \code{target} is input as a risk difference, and 
  transformed internally.  When the \code{fit} is of class \code{emaxalt} or
  \code{emaxsimobj}, the \code{binary} status is taken from the object and
  \code{binary} is ignored.}
}

\value{

Returns a vector with two elements:
  \item{targetDose}{The MLE of the dose achieving the target.}
  \item{seTD}{SE for target.dose}
}
\note{
Asymptotic SE computed using the delta method
}

\author{Neal Thomas }

\seealso{\code{\link{SeEmax}}, \code{\link{emaxalt}} }
\examples{
    \dontrun{

		## emaxsim changes the random number seed
    doselev<-c(0,5,25,50,100,250)
    n<-c(78,81,81,81,77,80)
    dose<-rep(doselev,n)

    ### population parameters for simulation
    e0<-2.465375 
    ed50<-67.481113 
    
   	dtarget<-100
		diftarget<-9.032497
		emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

    sdy<-7.967897
    pop<-c(led50=log(ed50),emax=emax,e0=e0)    
    meanresp<-emaxfun(dose,pop)  
    y<-rnorm(sum(n),meanresp,sdy)
    nls.fit<-nls(y ~ e0 + (emax * dose)/(dose + exp(led50)), 
                              start = pop, control = nls.control(
                              maxiter = 100),trace=TRUE,na.action=na.omit)

    targetD(nls.fit,10,modType=3)

    ###
    ### apply targetD to an emaxsim object
    ###
    nsim<-50
    sdy<-25
    gen.parm<-FixedMean(n,doselev,emaxfun(doselev,pop),sdy)  
    D4 <- emaxsim(nsim,gen.parm,modType=4)
    summary(D4,testalph=0.05)

    out<-NULL
    for(i in 1:nsim){
        out<-rbind(out,targetD(D4[i],target=4))
    }
}
    \dontshow{

		## emaxsim changes the random number seed
    doselev<-c(0,5,25,50,100,250)
    n<-c(78,81,81,81,77,80)
    dose<-rep(doselev,n)

    ### population parameters for simulation
    e0<-2.465375 
    ed50<-67.481113 
    
 		dtarget<-100
		diftarget<-9.032497
		emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

    sdy<-7.967897
    pop<-c(led50=log(ed50),emax=emax,e0=e0)    
    meanresp<-emaxfun(dose,pop)  
    y<-rnorm(sum(n),meanresp,sdy)
    nls.fit<-nls(y ~ e0 + (emax * dose)/(dose + exp(led50)), 
                              start = pop, control = nls.control(
                              maxiter = 100),trace=TRUE,na.action=na.omit)

    targetD(nls.fit,10,modType=3)
}
}
\keyword{nonlinear}

