\name{plot.plotB}
\alias{plot.plotB}

\title{Plot Bayes dose response curve and dose group means  }
\description{Plot a dose response curve fit by Bayes MCMC methods
(with optional posterior interval bars).  Also plot dose
group means (with optional CI bars)
  
}
\usage{
\method{plot}{plotB}( x,
	plotDif= FALSE, plotMed= FALSE, 
	plotResid=FALSE, predict= TRUE,
	logScale=FALSE,
	xlim,
	xat=NULL,
	ylim,  
	xlab, 
	ylab,labac='Act Comp',shapeac=8,colac='red',  
	symbolLabel='Group',symbolShape=8,
	symbolColor='red',symbolSize=4, ...)
}

\arguments{
  \item{x}{\code{\link{plotB}} object output from function \code{plotB}.}
  \item{plotDif}{Plot difference between doses and placebo. 
  It is assumed the lowest dose is placebo. If \code{activeControl}, the
  difference is with the active control mean, and the active 
  controls are not plotted.}
  \item{plotMed}{If TRUE, model-based curves are medians rather than means.}
  \item{plotResid}{If TRUE, a plot of the residuals formed from the
  dose group means minus the posterior dose group means.}
  \item{predict}{When predict=TRUE, predictive intervals for sample dose group
  proportions are plotted.  They are gray-shaded bars.}
  \item{logScale}{If \code{TRUE}, log scale is used for dose.}
  \item{xlim}{x-axis limits }
  \item{xat}{The points at which tick-marks are to be drawn. Errors occur if the 
  points are outside the range of xlim. By default (when NULL) tickmark locations are    computed.}
  \item{ylim}{y-axis limits }
  \item{xlab}{x-axis label }
  \item{ylab}{y-axis label }
  \item{labac}{x-axis label for the active control group.}
  \item{shapeac}{Shape of the symbol for the active control group.}
  \item{colac}{Color of the symbol for the active control group.}
  \item{symbolLabel}{Label given to symbol in plot legend.}
  \item{symbolShape}{A character vector with names giving the shapes assigned
  to different levels of variable \code{symbol}.  If a single shape is 
  specified, it is replicated for all dose groups.  See package
  \code{ggplot2} for symbol mappings.}
  \item{symbolColor}{A character vector with names giving the colors assigned
  to different levels of variable \code{symbol}.  If a single color is 
  specified, it is replicated for all dose groups.  See package
  \code{ggplot2} for color mappings.}
  \item{symbolSize}{The size of the symbol for the dose group sample means.   
   Set \code{symbolSize=0} to supress plotting.}
  \item{\ldots}{Additional parameters (not used)}
}
\value{ggplot object of the dose response curve, which will be plotted by
	default unless the output of the plot is assigned.  When an active
	control group is present, the value returned is an invisible list with the
	ggplot for the dosing data, and a second ggplot for the ac data.}
\details{
Produce additional plots from output of \code{\link{plotB}} without
any re-computing.  A plot is produced by default on return from the
function.  When active control is specified, the plot is 'printed'
within the function.  If there is a symbol group variable, it must be
specified when \code{plotB} is executed.  The symbol label, shape, color,
and size must be re-specified in subsequent plot requests.}

\note{
PlotB can also be used with draws from a prior distribution to
evaluate the prior dose response curve.
}

\author{ Neal Thomas }

\seealso{\code{\link{plotB}}, \code{\link{plotD}},  \code{\link{plot.fitEmax}} }
\examples{
\dontrun{
data("metaData")
exdat<-metaData[metaData$taid==6 & metaData$poptype==1,]

prior<-emaxPrior.control(epmu=0,epsca=100,difTargetmu=0,difTargetsca=100,dTarget=80.0,
        p50=3.75,sigmalow=0.01,sigmaup=20)
mcmc<-mcmc.control(chains=3)

msSat<-sum((exdat$sampsize-1)*(exdat$sd)^2)/(sum(exdat$sampsize)-length(exdat$sampsize))
fitout<-fitEmaxB(exdat$rslt,exdat$dose,prior,modType=4,
				count=exdat$sampsize,msSat=msSat,mcmc=mcmc)
parms<-coef(fitout)[,1:4]  #use first intercept

outB<-plotB(exdat$rslt,exdat$dose,parms, sigma2=(sigma(fitout))^2,
			ylab="Change in EDD")
			
plot(outB,plotDif=TRUE)
}
\dontshow{
data("metaData")
exdat<-metaData[metaData$taid==6 & metaData$poptype==1,]

parms<-matrix(c(1.611445,0.5160108,-77.73158,8.981834,
2.606269,0.3915947, -94.49052, 7.518533,
1.312385,0.9194756, -68.27381, 9.848746,
2.193071,0.5296447, -81.90235, 2.330520,
3.798612,0.3540451,-129.37529, 8.739098),ncol=4,byrow=TRUE)
sig2<-c(78.87238,87.05488,94.88178,78.54260,89.51093,79.74265)

outB<-plotB(exdat$rslt,exdat$dose,parms, sigma2=sig2,
			ylab="Change in EDD")
plot(outB,plotDif=TRUE)
}
}
\keyword{nonlinear}
\keyword{hplot}
