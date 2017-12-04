\name{FixedMean}
\alias{FixedMean}

\title{ Fixed means (proportions) random data constructor for emaxsim
 for continuous or binary data }
\description{
Creates a list object that contains inputs and  a function to create
simulated data sets with a common mean (proportion) for use in emaxsim
with normal or continuous data
}
\usage{
FixedMean(n, doselev, meanlev, resSD, parm = NULL, binary=FALSE)
}

\arguments{
  \item{n}{ Sample size for each dose group }
  \item{doselev}{Dose levels (including 0 for placebo) in the
study corresponding to \code{n}.  Must be in increasing order.}
  \item{meanlev}{ Mean response at each doselev }
  \item{resSD}{ Standard deviation for residuals within each dose
group (assumed common to all dose groups) }
  \item{parm}{ Population parameters that are 
saved for later reference, but are not used when creating simulated
   data.  \code{parm} can contain
parameters for a 3- or 4- parameter Emax model that generated
\code{meanlev}.  They should be stored in the order given in
    \code{\link{emaxfun}}.  Default is \code{NULL}. }
  \item{binary}{Normal data with homogeneous variance are generated unless
	  \code{binary} is \code{TRUE}, and then means are interpreted as
	  proportions and \code{0/1} data are generated.}
}

\value{
A list of length \code{2}.  
The first element is itself a list named \code{genP} that contains named elments
\code{n}, \code{resSD}, \code{doselev}, \code{dose}, \code{parm}, 
\code{binary}, and the
element \code{meanlev}, which is specific to \code{FixedMean}.  The second
element is a function named \code{genFun} that takes
\code{genP} as input and returns a list with named elements \code{meanlev}, 
\code{parm}, \code{resSD}, \code{y}.
}

\author{ Neal Thomas }

\seealso{ \code{\link{emaxsim}}, \code{\link{RandEmax}} }
\examples{

\dontrun{
##  example changes the random number seed

doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)  
  
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy,pop)  

### 4-parameter example

est<-c( log(6.67636), 2,0, -3.18230)
doselev<-c(0,5,10,25,50,150)
meanlev<-emaxfun(doselev,est)

gen.parm4<-FixedMean(n=c(99,95,98,94,98,99),doselev,
                     meanlev,resSD=3.87,parm=est)
D4 <- emaxsim(nsim=100,gen.parm4,modType=4,negEmax=TRUE)
summary(D4)
}
\dontshow{
##  example changes the random number seed

doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)  
  
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy,pop)  

}
}
\keyword{nonlinear}
