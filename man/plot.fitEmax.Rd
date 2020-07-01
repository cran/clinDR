\name{plot.fitEmax}
\alias{plot.fitEmax}
\title{Plot a Emax model and dose group means.
}
\description{
Plot an Emax model stored in an object
created by function \code{fitEmax}.
}
\usage{
\method{plot}{fitEmax}(
          x,int=0,plotResid=FALSE,clev=0.9,
          predict=TRUE,plotci=TRUE,plotDif=FALSE,
          xlab='Dose',
          ylab=ifelse(plotResid,'Residuals',ifelse(plotDif,
                'Difference With Placebo','Response')),
          symbol=NULL,symbolLabel='Group',symbolShape=8,
          symbolColor='red',symbolSize=4,
          bwidth=NULL,
          xlim=NULL,
          xat=NULL,
          ylim=NULL,
          logScale=FALSE,       
          ngrid=200,
          plot=TRUE, ...)
}

\arguments{
  \item{x}{Output of \code{\link{fitEmax}} with class
"fitEmax".}
  \item{int}{The index for the protocol (intercept) to use for the predictions 
  and computation of dose group means and standard errors.  The default
  value is \code{0}, which displays all protocols in a grid layout.}
  \item{plotResid}{If \code{TRUE}, a residual plot of the observed dose group means  
  is produced instead of a dose response curve plot.}
  \item{clev}{Confidence level for intervals about the estimated
mean for each dose.}
  \item{predict}{When predict=TRUE, predictive intervals for sample dose group
  means are plotted.  They are gray-shaded bars. If there is \code{>1} symbol group
  mean for a protocol/dose combination, then the smaller sample size is used when computing
  the prediction interval.}
  \item{plotci}{When plotCI=TRUE, confidence intervals for the population dose group
  means are plotted.  They are black bars.}
  \item{plotDif}{Plot difference between doses and placebo.  It is assumed the lowest
  dose in each protocol is placebo.}
  \item{xlab}{Label for the x-axis}
  \item{ylab}{Label for the y-axis}
  \item{symbol}{An optional grouping variable.  
  The values of symbol must correspond to the original data used in 
  \code{fitEmax}.}
  \item{symbolLabel}{Label given to symbol in plot legend.}
  \item{symbolShape}{A character vector with named elements giving 
  the shapes assigned
  to different levels of variable \code{symbol}.  If a single shape is 
  specified, it is replicated for all dose group means.  See package
  \code{ggplot2} for symbol mappings.}
  \item{symbolColor}{A character vector with named elements 
  giving the colors assigned
  to different levels of variable \code{symbol}.  If a single color is 
  specified, it is replicated for all dose group means.  See package
  \code{ggplot2} for color mappings.}
  \item{symbolSize}{The size of the symbol for the dose group sample means.   
   Set \code{symbolSize=0} to supress plotting the means.}
  \item{bwidth}{Width of the cap on the predictive interval bars.}
  \item{xlim}{Plot limits for the x-axis}
  \item{xat}{The points at which tick-marks are to be drawn. Errors occur if the 
  points are outside the range of xlim. By default (when NULL) tickmark locations are    computed.}  
  \item{ylim}{Plot limits for the y-axis}
  \item{logScale}{If \code{TRUE}, log scale is used for dose.}
  \item{ngrid}{The number doses evaluated when plotting the curve.}
  \item{plot}{Return plotting output without plotting.}
  \item{...}{No additional plotting options are currently used.}
}

\value{
A list with ggplot object, and a matrix with the confidence and prediction interval
limits.
}

\details{
Model estimates, standard errors, and confidence bounds are computed 
using function \code{\link{SeEmax}}.   

The function generates random numbers when \code{predict=TRUE}, so the random number generator/seed must
be set before the function is called for exact reproducibility.
}
\author{Neal Thomas }

\seealso{\code{\link{nls}}}
\examples{
### example changes the random number seed

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 

dtarget<-100
diftarget<-9.032497
emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

sdy<-8.0
pop.parm<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop.parm)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)

plot(testout)
}

\keyword{nonlinear}
