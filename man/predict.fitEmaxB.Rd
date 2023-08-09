\name{predict.fitEmaxB}
\alias{predict.fitEmaxB}
\title{
Estimated mean and posterior intervals
derived from a Bayesian hyperbolic or sigmiodial Emax model.
}
\description{
The mean/proportion response for different doses estimated from a Bayesian Emax model is
computed along with corresponding posterior intervals.  The results are computed
for a vector of input dose levels. The estimates are posterior means or medians of the
MCMC generated means/proportions.  
For binary outcomes, the estimated response rates are
computed on the logit scale and then back-transformed before forming the 
estimates and posterior intervals.
}
\usage{
\method{predict}{fitEmaxB}(object, dosevec, clev = 0.9,
	int = 1, dref = 0, xvec=NULL, ...)
}
\arguments{
  \item{object}{Output of \code{\link{fitEmax}} with class
"fitEmaxB".}
  \item{dosevec}{Vector of doses to be evaluated.}
  \item{clev}{Level for the posterior intervals about the 
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
Results computed from simple tabulations of the MCMC parameters
evaluated in the Emax function.

If baseline covariates were included in the fit and \code{xvec} is not specified, 
then the predicted value is the mean of the predictions for all patients in the
specified protocol.   Note that the protocol must be specified, or the 
prediction defaults to patients from the first protocol.  Note that for binary
data, the distinction between the mean of the predicted values and the predicted
value at the mean of the covariates can be important.
}
\value{
A list with estimated mean/proportion \code{(pred, predMed)}, lower bound, upper
bound, posterior SD, and corresponding values for differences
with the reference dose.  One value for each dose in dosevec.
The MCMC response means (proportions for binary data) are in
\code{simResp}, and the residual SD for continuous data are in
\code{sigsim}.}

\author{Neal Thomas }


\seealso{
\code{fitEmaxB}
}
\examples{

\dontrun{
data("metaData")
exdat<-metaData[metaData$taid==6 & metaData$poptype==1,]

prior<-emaxPrior.control(epmu=0,epsca=10,difTargetmu=0,difTargetsca=10,dTarget=80.0,
        p50=3.75,sigmalow=0.01,sigmaup=20)
mcmc<-mcmc.control(chains=3)

msSat<-sum((exdat$sampsize-1)*(exdat$sd)^2)/(sum(exdat$sampsize)-length(exdat$sampsize))
fitout<-fitEmaxB(exdat$rslt,exdat$dose,prior,modType=4,
				count=exdat$sampsize,msSat=msSat,mcmc=mcmc)

predout<-predict(fitout,dosevec=sort(unique(exdat$dose)))
}
}

