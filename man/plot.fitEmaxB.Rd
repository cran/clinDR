\name{plot.fitEmaxB}
\alias{plot.fitEmaxB}
\title{Plot a Emax model and dose group means.
}
\description{
Plot an Emax model stored in an object
created by function \code{fitEmaxB}.
}
\usage{
\method{plot}{fitEmaxB}(
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
  \item{x}{Output of \code{\link{fitEmaxB}} with class
"fitEmaxB".}
  \item{int}{The index for the protocol (intercept) to use for the predictions 
  and computation of dose group means/proportions.  The default
  value is \code{0}, which displays all protocols in a grid layout.}
  \item{plotResid}{If \code{TRUE}, a residual plot of the observed 
  dose group means/proportions  
  less the model-based MCMC median estimates of the means/proportions.}
  \item{clev}{Level for posterior probability intervals about the 
mean/proportion for each dose.}
  \item{predict}{When predict=TRUE, predictive intervals for sample dose group
  means/proportions are plotted.  They are gray-shaded bars. If there is \code{>1} symbol group
  mean/proportion for a protocol/dose combination, then the smaller sample size is used when computing
  the prediction interval.}
  \item{plotci}{When plotCI=TRUE, posterior intervals for the population dose group
  means/proportions are plotted.  They are black bars.}
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
  specified, it is replicated for all dose group means/proportions.  See package
  \code{ggplot2} for symbol mappings.}
  \item{symbolColor}{A character vector with named elements 
  giving the colors assigned
  to different levels of variable \code{symbol}.  If a single color is 
  specified, it is replicated for all dose group means/proportions.  See package
  \code{ggplot2} for color mappings.}
  \item{symbolSize}{The size of the symbol for the dose group sample means.   
   Set \code{symbolSize=0} to supress plotting the means.}
  \item{bwidth}{Width of the cap on the predictive interval bars.}
  \item{xlim}{Plot limits for the x-axis}
  \item{xat}{The points at which tick-marks are to be drawn. Errors occur if the 
  points are outside the range of xlim. By default (when NULL) tickmark locations are computed.}  
  \item{ylim}{Plot limits for the y-axis}
  \item{logScale}{If \code{TRUE}, log scale is used for dose.}
  \item{ngrid}{The number doses evaluated when plotting the curve.}
  \item{plot}{Return plotting output without plotting.}
  \item{...}{No additional plotting options are currently used.}
}

\value{
A list with ggplot object, and posterior and prediction interval
limits.
}

\details{
Model-based medians, standard deviations, and interval bounds for the dose groups means/proportions based on the MCMC parameters evaluated in the Emax function.   

The function generates random numbers when \code{predict=TRUE}, so the random number generator/seed must
be set before the function is called for exact reproducibility.

If baseline covaraties were included in the fit, then the mean of the predictions for
the protocol given by \code{int} is plotted.  This can be computationally intensive
when the dosing grid is dense, the MCMC sample size is large, and the input sample
size is large.  Consider reducing \code{ngrid} in this situation.  
Note that the protocol must be specified, or the prediction defaults to patients
from the first protocol.
}
\author{Neal Thomas }

\seealso{\code{\link{fitEmaxB}}}
\examples{
\dontrun{

data("metaData")
exdat<-metaData[metaData$taid==1,]

prior<-emaxPrior.control(epmu=0,epsca=4,difTargetmu=0,difTargetsca=4,dTarget=20,
				p50=(2+5)/2,
				sigmalow=0.01,sigmaup=3)
										
mcmc<-mcmc.control(chains=3)

msSat<-sum((exdat$sampsize-1)*(exdat$sd)^2)/(sum(exdat$sampsize)-length(exdat$sampsize))
fitout<-fitEmaxB(exdat$rslt,exdat$dose,prior,modType=4,prot=exdat$protid,
				count=exdat$sampsize,msSat=msSat,mcmc=mcmc)
plot(fitout)
}
}

\keyword{nonlinear}
