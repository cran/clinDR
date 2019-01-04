\name{RandEmax}
\alias{RandEmax}

\title{Random data constructor function for emaxsim creating random parameters
for an Emax model for continuous or binary data.}
\description{
Creates a list object that contains inputs and a function to create
simulated data sets for emaxsim.  Data sets are created by
generating random parameters from beta or log-normal distributions for 
a 3/4 parameter Emax model. For binary data, the Emax model is on the logit scale
and then back-transformed.
}
\usage{
RandEmax(n, doselev, 
	parmEmax,
	parmE0,
	p50,
	parmED50=c(3,0.79,0.6),
	parmLambda=c(3.03,18.15,0,6),
	resSD,
	dfSD=Inf,
	binary=FALSE)
}
\arguments{
  \item{n}{Sample size for each dose group. }
  \item{doselev}{Dose levels (including 0 for placebo) included in the
study corresponding to \code{n}.  Must be in increasing order.}
  \item{parmEmax}{Vector with mean and standard deviation for a random
normal Emax}  
  \item{parmE0}{Vector with mean and standard deviation for a random
normal intercept.}
  \item{p50}{The predicted ED50}
  \item{parmED50}{  The log(ED50) is generated from a t-distribution
  with \code{df=parmED50[1]}, \code{mean=log(p50)+parmED50[2]}, 
  and \code{scale=parmED50[3]}.
  The default values are taken from the reference below.
  }
  \item{parmLambda}{For a beta distributed sigmoid lambda, a vector with
(df1,df2,lower bound, upper bound).  For a hyperbolic model, lambda=1.  }
  \item{resSD}{Standard deviation for residuals within each dose (normal data only)}
  \item{dfSD}{If a finite value is specified, the within-dose group SD is randomly
  generated from resSD times sqrt(dfSD/chisquare(dfSD))), which is the form of a posterior
  distribution for a SD based on a existing sample.}
  \item{binary}{When \code{TRUE}, \code{0/1} data are generated from the
	  Emax model, which is computed on the logit scale and then 
	  backtransformed to yield proportions.}
}
\details{
All parameters are independent.  Normal data are
generated from the dose response curves with homogeneous-variance
normal residuals.  Binary data are 0/1 generated from Bernoulli 
distributions with proportions computed by transforming the Emax
model output from the logit to proportion scale.  
Default values are based on recommendations in

Thomas, N., Sweeney, K., and Somayaji, V. (2014).  Meta-analysis 
of clinical dose response in a large drug development portfolio.  
<doi:10.1080/19466315.2014.924876>
}
\value{
A list of length \code{2}.  
The first element is itself a list named \code{genP} that contains named elments
\code{n}, \code{resSD}, \code{dfSD}, \code{doselev}, \code{dose},  
\code{binary} and the
elements \code{parmE0}, \code{p50}, \code{parmED50}, \code{parmEmax}, 
and \code{parmLambda}. 
which are specific to \code{RandEmax}.  The second
element is a function named \code{genFun} that takes
\code{genP} as input and returns a list with named elements \code{meanlev}, 
\code{parm}, \code{resSD}, \code{y}.  
}

\author{ Neal Thomas }

\seealso{ \code{\link{emaxsim}}, \code{\link{FixedMean}} }

\examples{
simParm<-RandEmax(n=c(99,95,98,94,98,98),doselev=c(0,5,10,25,50,150),
         parmE0=c(-2.6,2.5),p50=25,parmEmax=c(-1.25,2),resSD=3.88)
}
\keyword{nonlinear}


