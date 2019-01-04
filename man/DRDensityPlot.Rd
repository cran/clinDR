\name{DRDensityPlot}
\alias{DRDensityPlot}

\title{Plot Bayes or confidence interval density contours over a grid of points (usually dose or time)} 
\description{Density plot for distributions conditional on a variable. 
A grid of values are specified for the conditioning variable, which
is plotted on the horizontal axis.  The conditioning variable is 
typically dose or time}
  
\usage{
DRDensityPlot(x,qL,qH,qlevL=c(0.025,0.05,0.10,0.25),
		xlim,ylim,xlab='x',ylab='y')
}

\arguments{
  \item{x}{A grid of conditioning values to be plotted on the horizontal axis.  This grid typically
  represents dose or time.}
  \item{qL}{Lower percentiles, confidence or probabiity levels.  \code{qL} is a matrix with rows 
  corresponding to \code{x}, and columns corresponding to \code{qlevL}.  The percentiles must be 
  increasing in order and less that 0.50.}
  \item{qH}{Upper percentiles, confidence or probabiity levels.  \code{qH} levels correspond to
  the \code{qL} levels but are ordered from highest to lowest (1-\code{qlevL}), with the 
  smallest greater than 0.50.}
  \item{qlevL}{Density intervals are formed with percentile boundaries at (qlevL,1-qlevL).  
  \code{qlevL} must be increasing between (0,0.5).}
  \item{xlim}{Plot limits for the x-axis}
  \item{ylim}{Plot limits for the y-axis}
  \item{xlab}{x-axis label }
  \item{ylab}{y-axis label }
 
}
\details{
The function takes as input percentiles defining confidence 
intervals or Bayesian probability 
intervals at different levels (e.g. 5,95, 25,75) for distributions conditional
on a variable that is typically dose or time.  Regions defined by different
confidence/probability levels are represented by different levels of shading.
The input parameter, \code{qlevL}, is used only to define the input in the
matrices \code{qL} and \code{qH}.  The \code{qlevL} is not used for any numerical 
calculations, which must be done before executing the function.
}
\value{
Plotted output only.
}

\author{ Neal Thomas }

\seealso{\code{\link{plotBdensity}}}
\examples{

\dontrun{
data(examples14)
exdat<-examples14[[32]]

fitout<-fitEmax(exdat$y,exdat$dose,modType=3,count=exdat$nsize,
								msSat=(exdat$sd)^2)

dgrid<-seq(0,1,length=100)
seout95<-predict(fitout,dgrid,clev=0.95)
seout90<-predict(fitout,dgrid,clev=0.9)
seout80<-predict(fitout,dgrid,clev=0.8)
seout50<-predict(fitout,dgrid,clev=0.5)

qlev<-c(0.025,0.05,0.10,0.25)

qL<-cbind(seout95$ubdif,seout90$ubdif,seout80$ubdif,seout50$ubdif)
qH<-cbind(seout95$lbdif,seout90$lbdif,seout80$lbdif,seout50$lbdif)

DRDensityPlot(dgrid,qL,qH,qlevL=qlev,xlab='Dose',ylab='Diff with PBO')

}
\dontshow{
data(examples14)
exdat<-examples14[[32]]

dgrid<-seq(0,1,length=5)

qlev<-c(0.10,0.25)

qL<-matrix(c(0.000000, 0.000000,
1.181590, 1.093189,
1.301505, 1.220726,
1.354046, 1.273955,
1.384266, 1.303586),ncol=2,byrow=TRUE)

qH<-matrix(c(0.0000000, 0.0000000,
0.8083449, 0.8967468,
0.9604440, 1.0412232,
1.0158898, 1.0959808,
1.0436238, 1.1243036),ncol=2,byrow=TRUE)

DRDensityPlot(dgrid,qL,qH,qlevL=qlev,xlab='Dose',ylab='Diff with PBO')
}
}
\keyword{nonlinear}
\keyword{hplot}
