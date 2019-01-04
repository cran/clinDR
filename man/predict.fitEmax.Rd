\name{predict.fitEmax}
\alias{predict.fitEmax}
\title{Estimated mean/proportion and confidence intervals
derived from the maximum likelihood fit of a 3- or 4- parameter Emax model.
}
\description{
The estimated means from an Emax model is
computed along with confidence bounds.  The results are computed
for a vector of input dose levels.  For binary outcomes, the results are
computed on the logit scale and then back-transformed.
}
\usage{
\method{predict}{fitEmax}(object,dosevec,clev=0.9, 
                  int=1,dref=0, xvec=NULL, ...)
}

\arguments{
  \item{object}{Output of \code{\link{fitEmax}} with class
"fitEmax".}
  \item{dosevec}{Vector of doses to be evaluated.}
  \item{clev}{Confidence level for intervals about the estimated
mean/proportion at each dosevec.}
  \item{int}{The index for the protocol (intercept) to use for the predictions}
    \item{dref}{Differences in response between \code{doselev}
  and \code{dref} are computed.}
  \item{xvec}{The vector of centered baseline values for the prediction model when
  \code{xbase} was specified in the model fit.  Centering must be done using the
  protocol-specific means consistent with \code{int}.  See details for the default
  calculations when \code{xvec} is not specified.}  
  \item{...}{No additonal parameters will be utilized.}
}

\details{
Model estimates, standard errors, and confidence bounds are computed with the 
function \code{\link{SeEmax}}. 

If baseline covariates were included in the fit and \code{xvec} is not specified, 
then the predicted value is the mean of the predictions for all patients in the
specified protocol.   Note that the protocol must be specified, or the 
prediction defaults to patients from the first protocol.  Note that for binary
data, the distinction between the mean of the predicted values and the predicted
value as the mean of the covariates can be important.
}
\value{
A list with estimated dose group means/proportions, lower bound, upper
bound, SE, and corresponding values for differences
with the reference dose.  One value for each dose in dosevec.
}

\author{Neal Thomas }

\seealso{\code{\link{nls}}}
\examples{

\dontrun{
## this example changes the random number seed
doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop.parm<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop.parm)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)
predout<-predict(testout,dosevec=c(20,80),int=1)
}
\dontshow{
## this example changes the random number seed
doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop.parm<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop.parm)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)
predout<-predict(testout,dosevec=c(20,80),int=1)
}
}

\keyword{nonlinear}
