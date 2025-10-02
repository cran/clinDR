\name{bpchkMonoEmax}
\alias{bpchkMonoEmax}

\title{Bayes posterior predictive test for Emax (monotone) model fit}
\description{Bayes posterior predictive test for an Emax (monotone) model  
fit comparing the best response from lower doses to
the response from the highest dose.  \code{checkMonoEmax} is deprecated.  
See \code{bpchkMonoEmax}.}
  
\usage{
bpchkMonoEmax(x, trend='positive', protSel=1)
}

\arguments{
  \item{x}{Output object of class 'fitEmaxB'.}
  \item{trend}{The default is 'positive', so high values for lower doses
  yield small Bayesian predictive probabilities.  Set \code{trend} to 'negative'
  for dose response curves with negative trends.}
  \item{protSel}{The test is applied to the data from a single protocol.  
  The protocol can be selected if the model was fit to data from more 
  than one protocol.  The \code{protSel} must match a 
  protocol value input to \code{fitEmaxB} or it numerical index value, 1,2,...}
}
\details{
The Bayesian predictive p-value 
is the posterior probability that a dose group
sample mean in a new study with the same sample sizes would yield a 
higher (or lower for negative trend) difference for one of the lower doses versus 
the highest dose than was actually obtained from the real sample.  There must
be at least two non-placebo dose groups (NA returned otherwise).  Placebo
response is excluded from the comparisons.
  
The function generates random numbers, so the random number generator/seed must
be set before the function is called for exact reproducibility.

When \code{fitEmaxB} is applied to first-stage fitted model output with a
non-diagonal variance-covariance matrix, the predictive draws are selected
from a multivariate model with means computed from the MCMC-generated
parameters and input asymptotic variance-covariance matrix \code{vcest}.  
If the fitted model was applied to binary data, the GOF statistic is 
computed based on the logit rather than observed dose group sample 
proportion scale.  This differs from the setting with patient-level
data input to \code{fitEmaxB}.
}
\value{
Returns a scalar Bayesian predictive p-value.
}

\references{
Thomas, N., Sweeney, K., and Somayaji, V. (2014).  Meta-analysis of clinical dose response in a large drug development portfolio, Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317.
<doi:10.1080/19466315.2014.924876>

Thomas, N., and Roy, D. (2016).  Analysis of clinical dose-response in small-molecule drug development:  2009-2014. Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317
<doi:10.1080/19466315.2016.1256229>

Wu, J., Banerjee, A., Jin, B. Menon, M. S., Martin, S. and Heatherington, A. (2017). Clinical dose response for a broad set of biological products: A model-based meta-analysis. Statistical Methods in Medical Research. <doi:10.1177/0962280216684528>
}

\author{ Neal Thomas }

\seealso{\code{\link{plot.plotB}}, \code{\link{plotD}},  \code{\link{plot.fitEmax}} }

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
parms<-coef(fitout)[,1:4]  #use first intercept

checkMonoEmax(fitout, trend='negative')
      
}
}
\keyword{nonlinear}

