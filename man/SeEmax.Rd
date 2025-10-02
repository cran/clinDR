\name{SeEmax}
\alias{SeEmax}

\title{ Asymptotic SE for dose response estimates from a
3- or 4- parameter Emax model }
\description{
  Compute the asymptotic SE for dose response estimates based on the 
asymptotic variance-covariance matrix from the fit of a 3- or 4-parameter
Emax model}
\usage{
SeEmax(fit, doselev, modType, dref=0, nbase=0, x=NULL, 
				binary=FALSE, clev=0.9)
}

\arguments{
  \item{fit}{Output of \code{\link{nls}} fit to a 3- or 4-parameter Emax
model.  The order of the parameters in the fit must be (log(ed50),emax,e0) or
(log(ed50),lambda,emax,e0).
 Alternatively, fit can be a list with the first element the coefficient
vector, and the second element the variance-covariance matrix.  List input
can be used with multiple protocols and baseline covariates (see details).}
  \item{doselev}{SEs are evaluated at vector of doses }
  \item{modType}{modType=3,4 for a 3 or 4 parameter model.}
  \item{dref}{A reference dose (0 by default) for contrasts, but other values 
   can be specified.  If specified, a single reference value must be given.
   }
   
  \item{nbase}{The number of baseline predictors included in the model.}
  \item{x}{The model is evaluated at baseline covariate values, \code{x}.
  If \code{x} is a matrix, then each row is a vector of baseline 
  predictors, and the results are for the dose response averaged over
  all of the predictors in \code{x}.}
  \item{binary}{Emax model on logistic scale, then backtransformed.}
  \item{clev}{Confidence level for intervals.}
}

\value{
Returns a list:
  \item{doselev}{Doses to evaluate}
  \item{dref}{Differences in response between \code{doselev}
  and \code{dref} are computed.}
  \item{fitpred }{Estimated dose response at doselev}
  \item{sepred}{SE for estimated dose responses}
  \item{fitdif}{Estimated response at doselev minus estimated
       response at placebo}
   \item{sedif}{SE for fitdif estimated differences}
   \item{fitref}{Estimated dose response at the reference dose.}
   \item{seref}{SE for the estimated dose response at the reference
   dose}
   \item{covref}{The covariance between each estimated response and
   the estimated response at the reference dose.  These covariances can
   be used to compute asymptotic variances of differences after 
   back-transformation (e.g., for logistic regression with binary data).
   }
}
\details{
The Emax models supported by \code{SeEmax} should now be fit using 
\code{fitEmax} and \code{predict.fitEmax}.  \code{SeEmax} remains available
primarily for backward compatibility.

\code{SeEmax} can be used with models that allow different placebo response
for multiple protocols by selecting the intercept for a specific protocol.  
Coeficients for baseline covariates can also be included following the intercept.
The variance-covariance matrix from the full model must be subsetted to match
the included coeficients (i.e., the rows and columns corresponding to the
omitted intercepts must be removed).  List input must be used for the more
general models. 
}
\references{
Bates, D. M. and Watts, D. G. (1988) Nonlinear Regression
     Analysis and Its Applications, Wiley
 }
\author{Neal Thomas }
\seealso{\code{fitEmax}}
\examples{

\dontrun{

## this example changes the random number seed
doselev<-c(0,5,25,50,100,250)
n<-c(78,81,81,81,77,80)
dose<-rep(doselev,n)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
led50<-log(ed50)
lambda=1.8

dtarget<-100
diftarget<-9.032497
emax<-solveEmax(diftarget,dtarget,log(ed50),lambda,e0)

sdy<-7.967897
pop<-c(led50=led50,lambda=lambda,emax=emax,e0=e0)    
meanresp<-emaxfun(dose,pop)  
y<-rnorm(sum(n),meanresp,sdy)
nls.fit<-nls(y ~ e0 + (emax * dose^lambda)/(dose^lambda + exp(led50*lambda)), 
                         start = pop, control = nls.control(
                         maxiter = 100),trace=TRUE,na.action=na.omit)


SeEmax(nls.fit,doselev=c(60,120),modType=4)
SeEmax(list(coef(nls.fit),vcov(nls.fit)),c(60,120),modType=4)
}
\dontshow{

## this example changes the random number seed
doselev<-c(0,5,25,50,100,250)
n<-c(78,81,81,81,77,80)
dose<-rep(doselev,n)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
led50<-log(ed50)

dtarget<-100
diftarget<-9.032497
lambda=1.8
emax<-solveEmax(diftarget,dtarget,log(ed50),lambda,e0)

sdy<-7.967897
pop<-c(led50=led50,lambda=lambda,emax=emax,e0=e0)    
meanresp<-emaxfun(dose,pop)  
y<-rnorm(sum(n),meanresp,sdy)
nls.fit<-nls(y ~ e0 + (emax * dose^lambda)/(dose^lambda + exp(led50*lambda)), 
                         start = pop, control = nls.control(
                         maxiter = 100),trace=TRUE,na.action=na.omit)


SeEmax(nls.fit,doselev=c(60,120),modType=4)
}
}

\keyword{nonlinear}

