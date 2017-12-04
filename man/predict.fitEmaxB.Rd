\name{predict.fitEmaxB}
\alias{predict.fitEmaxB}
\title{
Estimated mean and posterior intervals
derived from a Bayesian hyperbolic or sigmiodial Emax model.
}
\description{
The mean response for different doses estimated from a Bayesian Emax model is
computed along with corresponding posterior intervals.  The results are computed
for a vector of input dose levels. The estimates are posterior means of the
dose response means.  
For binary outcomes, the estimated response rates are
computed on the logit scale and then back-transformed before forming the 
estimates and posterior intervals.
}
\usage{
\method{predict}{fitEmaxB}(object, dosevec, clev = 0.9, 
								int = 1, dref = 0, ...)
}
\arguments{
  \item{object}{Output of \code{\link{fitEmax}} with class
"fitEmaxB".}
  \item{dosevec}{Vector of doses to be evaluated.}
  \item{clev}{Posterior intervals about the 
proportion of responders at each dosevec.}
  \item{int}{The index for the protocol (intercept) to use for the predictions}
    \item{dref}{Differences in response between \code{doselev}
  and \code{dref} are computed.}
  \item{...}{No additonal parameters will be utilized.}
}
\details{
Results computed from simple tabulations of the MCMC parameters
evaluated in the Emax function.
}
\value{
A list with estimated mean, lower bound, upper
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
data("examples14")
exdat<-examples14[[6]]

prior<-prior.control(epmu=0,epsd=10,emaxmu=0,emaxsd=10,p50=0.25,
				sigmalow=0.01,sigmaup=3)
mcmc<-mcmc.control(chains=3)

fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
				count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)

predout<-predict(fitout,dosevec=sort(unique(exdat$dose)))
}
}

