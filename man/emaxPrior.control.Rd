\name{emaxPrior.control}
\alias{emaxPrior.control}
\title{
Set the parameters of the prior distribution for the Emax model implemented in \code{fitEmaxB}.
}
\description{
Set the parameters of the prior distribution for the Emax model implemented in \code{fitEmaxB.}.
}
\usage{
emaxPrior.control(epmu=NULL,epsca=NULL,
	difTargetmu=NULL,difTargetsca=NULL,
	dTarget=NULL,p50=NULL,
	sigmalow=NULL,sigmaup=NULL,
	effDF=parmDF,parmDF=5,
	loged50mu=0.0,loged50sca=1.73,
	loglammu=0.0,loglamsca=0.425,parmCor=-0.45,
	basemu=NULL,basevar=NULL,binary=FALSE)
}
\arguments{
  \item{epmu}{
Mean for \code{E0} in a t-prior distribution.  Logistic scale for binary data.
}
  \item{epsca}{
The scale parameter for \code{E0} in a t-prior distribution.  Logistic scale for binary data.
}
  \item{difTargetmu}{
Mean for the prior distribution of the effect at dose \code{dTarget} versus placebo.
Logistic scale for binary data. 
}
  \item{difTargetsca}{
The scale parameter for the prior distribution of the effect at dose \code{dTarget} versus placebo.
Logistic scale for binary data.  
}
	\item{dTarget}{Target dose for prior effect.  Typically the highest dose planned
	and/or the proof-of-concept dose.}
  \item{p50}{
Projected \code{ED50}.  See references for its use in creating the prior distribution for the \code{ED50}.
}
  \item{sigmalow}{
Lower bound for a uniform prior distribution for the residual SD (continuous data).
}
  \item{sigmaup}{
Upper bound for a uniform prior distribution for the residual SD (continuous data).  
}
	\item{effDF}{The degrees of freedom for the log-t prior distributions for the \code{placebo}
	and \code{difTarget} parameters.  If a vector of length 2 is specified, the first value is the
	degrees of freedom for \code{placebo} and the second for \code{difTarget}.}
  \item{parmDF}{
The degrees of freedom of the bivariate log-t prior distribution for the \code{ED50} and \code{lambda} parameters.  
}
  \item{loged50mu}{
Mean of prior t-distribution for the \code{log(ED50)}.  See references for 
its default value and interpretation.
}
  \item{loged50sca}{
Scale (analogous to SD) of the prior t-distribution for the \code{log(ED50)}.  
}
  \item{loglammu}{
Mean of prior t-distribution for the Hill parameter lambda.  See 
references for its default value and interpretation.
}
  \item{loglamsca}{
Scale (analogous to SD) of the prior t-distribution for the Hill 
parameter lambda.  
}

	\item{parmCor}{Correlation for the bivariate log-t prior distribution for the \code{ED50} and \code{lambda} parameters.}

	\item{basemu}{
	A vector of prior means for the covariate regression parameters.
	}
	\item{basevar}{
	The prior variance-covariance matrix for the covariate regression parameters.
	The covariate regression parameters are apriori independent of the other dose 
	response model parameters.
	}
  \item{binary}{
Set to \code{TRUE} for binary data applications.  Used to check for consistency in usage.
}
}

\value{
List of class \code{emaxPrior} of prior parameter values for use in \code{fitEmaxB}.
}

\details{
The prior distribution is based on meta-analyses of dose 
response described in the references.  
The E0 and difTarget parameters have independent t-distribution prior 
distributions.  For binary data, these
parameters are computed on the logistic scale.  The prior means and scales
of these parameters must be assigned compound-specific values.  
The predicted ED50 at the study design stage must must also be
specified as 'P50'.   For
continuous data, the prior distribution for the residual SD is uniform on a 
user-specifed scale.

The prior distribution of the log(ED50) has a t-distribution
centered at log(P50), with scale, degrees of freedom (parmDF), 
and offset to the P50, 
defaulting to values given in the references (these can be changed, but they
are difficult to interpret outside the context of the meta-analyses).  
If \code{modType=4}, the prior distribution for the Hill parameter
is also t-distribution with parmDF degrees of freedom and corParm
correlation with the log(ED50). 
}


\references{
Thomas, N., Sweeney, K., and Somayaji, V. (2014).  Meta-analysis of clinical dose response in a large drug development portfolio, Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317.
<doi:10.1080/19466315.2014.924876>

Thomas, N., and Roy, D. (2016).  Analysis of clinical dose-response in small-molecule drug development:  2009-2014. Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317
<doi:10.1080/19466315.2016.1256229>

Wu, J., Banerjee, A., Jin, B., Menon, S., Martin, S., and 
Heatherington, A. (2017).  Clinical dose-response for a broad set of 
biological products: A model-based meta-analysis. Vol. 9, 2694-2721.   
<doi:10.1177/0962280216684528?>

}
\author{
Neal Thomas
}

\seealso{
\code{fitEmaxB}
}
