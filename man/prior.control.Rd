\name{prior.control}
\alias{prior.control}
\title{
Set the parameters of the prior distribution for the Emax model implemented in \code{fitEmaxB}.
}
\description{
Set the parameters of the prior distribution for the Emax model implemented in \code{fitEmaxB}.  \code{prior.control} is deprecated.  See \code{emaxPrior.control}.}
\usage{
prior.control(epmu = NULL, epsd = NULL, emaxmu = NULL, 
	emaxsd = NULL, p50 = NULL, 
	sigmalow = NULL, sigmaup = NULL, 
	led50mu = 0.79, led50sca = 0.6, edDF = 3, 
	lama = 3.03, lamb = 18.15, lamsca = 6, 
	basemu=NULL,basevar=NULL,
	binary = FALSE)
}
\arguments{
  \item{epmu}{
Mean for \code{E0} in a normal prior distribution.  Logistic scale for binary data.
}
  \item{epsd}{
SD for \code{E0} in a normal prior distribution.  Logistic scale for binary data.
}
  \item{emaxmu}{
Mean for \code{Emax} in a normal prior distribution.  Logistic scale for binary data. 
}
  \item{emaxsd}{
SD for \code{Emax} in a normal prior distribution.  Logistic scale for binary data.  
}
  \item{p50}{
Projected \code{ED50}.  See reference for its use in creating the prior distribution for the \code{ED50}.
}
  \item{sigmalow}{
Lower bound for a uniform prior distribution for the residual SD (continuous data).
}
  \item{sigmaup}{
Upper bound for a uniform prior distribution for the residual SD (continuous data).  
}
  \item{led50mu}{
Mean of log-t prior distribution for the \code{ED50} before final scaling.  See reference for its interpretation in the prior distribution for the \code{ED50}.
}
  \item{led50sca}{
Scale (analogous to SD) of the log-t prior distribution for the \code{ED50}.  
}
  \item{edDF}{
The degrees of freedom of the log-t prior distribution for the \code{ED50}.  
}
  \item{lama}{
  Parameter in the re-scaled beta distribution for Hill slope parameter in the sigmoidal
  Emax model.  See reference for it use and empirical basis.
}
  \item{lamb}{
  Parameter in the re-scaled beta distribution for Hill slope parameter in the sigmoidal
  Emax model.
}
  \item{lamsca}{
The beta prior distribution for the Hill parameter is re-scaled to have support on (0,lamsca).
}
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
List of prior parameter values for use in \code{fitEmaxB}.
}

\details{
The prior distributions are based two meta-analyses of dose response described 
in the references.  Each parameter is independent in the prior distribution.  The
E0 and Emax parameters have normal prior distributions.  For binary data, these
parameters are computed on the logistic scale.  The predicted ED50 must be
specified as 'P50'.  The prior distribution of the log(ED50) has a t-distribution
centered at log(P50), with scale, degrees of freedom, and offset to the P50, 
defaulting to values given in the references (these can be changed, but they
are difficult to interpret outside the context of the meta-analyses).  
If \code{modType=4}, the prior distribution for the Hill parameter
is a beta distribution scaled to
(0,lamsca).  The default degrees of freedom were obtained from the
meta-analyses.   For
continuous data, the prior distribution for the residual SD is uniform on a 
user-specifed scale.
}


\references{
Thomas, N., Sweeney, K., and Somayaji, V. (2014).  Meta-analysis of clinical dose response in a large drug development portfolio, Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317.
<doi:10.1080/19466315.2014.924876>

Thomas, N., and Roy, D. (2016).  Analysis of clinical dose-response in small-molecule drug development:  2009-2014. Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317
<doi:10.1080/19466315.2016.1256229>
}
\author{
Neal Thomas
}

\seealso{
\code{fitEmaxB}
}


