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
genp<-FixedMean(n,doselev,meanlev,sdy,pop)  

### binary example
n<-rep(500,5)
doselev<-c(0,5,25,50,1000)
dose<-rep(doselev,n)

e0<- qlogis(0.2)
ed50<-20
emax<-qlogis(0.6)-qlogis(0.2)
lambda<-2
pop<-c(log(ed50),lambda,emax,e0)
meanlev<-plogis(emaxfun(doselev,pop))

genp<-FixedMean(n,doselev,meanlev,sdy,pop,binary=TRUE)  

tapply(genp$genFun(genp$genP)$y,dose,mean)
meanlev

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
