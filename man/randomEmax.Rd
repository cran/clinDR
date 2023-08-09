\name{randomEmax}
\alias{randomEmax}

\title{Random data constructor function for emaxsim(B) 
creating random parameters
for an Emax model for continuous or binary data.}
\description{
Creates a list object that contains inputs and a function to create
simulated data sets for emaxsim(B).  Data sets are created by
generating random parameters from an \code{emaxPrior.control()} object for 
a 3/4 parameter Emax model. For binary data, the Emax model is 
on the logit scale and then back-transformed. 
}
\usage{
randomEmax(x,n,doselev,modType=c('4','3'))
}
\arguments{
	\item{x}{Object of type \code{emaxPrior} created by function
	\code{emaxPrior.control}, that specifies a prior distribution
	for the Emax model parameters.}
  \item{n}{Sample size for each dose group. }
  \item{doselev}{Dose levels (including 0 for placebo) included in the
study corresponding to \code{n}.  Must be in increasing order.}
	\item{modType}{Specifies a 4-parameter sigmoidal Emax model, or
	a 3-parameter hyperbolic Emax model}
}
\details{
Normal data are
generated from the dose response curves with homogeneous-variance
normal residuals.  Binary data are 0/1 generated from Bernoulli 
distributions with proportions computed by transforming the Emax
model output from the logit to proportion scale.  
Default values are based on recommendations in the references.
}

\references{
Thomas, N., Sweeney, K., and Somayaji, V. (2014).  Meta-analysis of clinical dose response in a large drug development portfolio, Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317.
<doi:10.1080/19466315.2014.924876>

Thomas, N., and Roy, D. (2016).  Analysis of clinical dose-response in small-molecule drug development:  2009-2014. Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317
<doi:10.1080/19466315.2016.1256229>

Wu, J., Banerjee, A., Jin, B. Menon, M. S., Martin, S. and Heatherington, A. (2017). Clinical dose response for a broad set of biological products: A model-based meta-analysis. Statistical Methods in Medical Research. <doi:10.1177/0962280216684528>
}

\value{
A list of length \code{2}.  
The first element is itself a list named \code{genP} that contains named elments
\code{n}, \code{doselev}, \code{dose},  
\code{modType} and the \code{emaxPrior} object \code{x}.  
The second element is a function named \code{genFun} that takes
\code{genP} as input and returns a list with named elements \code{meanlev}, 
\code{parm}, \code{resSD}, \code{y}.  
}

\author{ Neal Thomas }

\seealso{ \code{\link{emaxsimB}}, 
\code{\link{emaxsim}}, \code{\link{FixedMean}} }

\examples{
prior<-emaxPrior.control(epmu=0,epsca=4,
	difTargetmu=0,difTargetsca=4,dTarget=20,
	p50=(2+5)/2,
	sigmalow=0.01,sigmaup=3)
						
simParm<-randomEmax(x=prior,n=c(99,95,98,94,98,98),
	doselev=c(0,5,10,25,50,150),modType="4")
		
# D1 <- emaxsimB(nsim=10,simParm,prior,nproc=1)
}
\keyword{nonlinear}


