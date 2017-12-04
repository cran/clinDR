\name{plot.plotB}
\alias{plot.plotB}

\title{Plot Bayes dose response curve and dose group means  }
\description{Plot a dose response curve fit by Bayes MCMC methods
(with optional posterior interval bars).  Also plot dose
group means (with optional CI bars)
  
}
\usage{
\method{plot}{plotB}( x, 
		plotDif= FALSE, plotMed= FALSE, plotResid=FALSE,
		predict= TRUE,
		log=FALSE,
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
  \item{plotDif}{If TRUE, plot difference between doses and placebo.}
  \item{plotMed}{If TRUE, model-based curves are medians rather than means.}
  \item{plotResid}{If TRUE, a plot of the residuals formed from the
  dose group means minus the posterior dose group means.}
  \item{predict}{When predict=TRUE, predictive intervals for sample dose group
  proportions are plotted.  They are gray-shaded bars.}
  \item{log}{If \code{TRUE}, log scale is used for dose.}
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
data("examples14")
exdat<-examples14[[6]]

prior<-prior.control(epmu=0,epsd=10,emaxmu=0,emaxsd=10,p50=0.25,
				sigmalow=0.01,sigmaup=3)
mcmc<-mcmc.control(chains=3)

fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
				count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)
parms<-as.matrix(fitout$estanfit)

outB<-plotB(exdat$y,exdat$dose,cbind(parms[,"led50"],parms[,"lambda"],
      parms[,"emax"],parms[,"e0[1]"]), sigma2=(parms[,'sigma'])^2,
			ylab="Change in EDD")
			
plot(outB,plotDif=TRUE)
}
\dontshow{
data("examples14")
exdat<-examples14[[6]]

parms<-matrix(c(-0.1665350, 0.3657811, -5.660137, 1.744753, 0.4050860,
-0.8463137, 0.3837361, -4.877676, 1.784098, 0.3943782,
-1.1811274, 0.3767222, -4.921861, 1.873861, 0.4266011,
 0.4729616, 0.3157714, -6.322768, 1.780517, 0.3646588,
 0.4255880, 0.3336959, -6.251558, 1.775438, 0.3657461),ncol=5,byrow=TRUE)

outB<-plotB(exdat$y,exdat$dose,parms[,1:4], sigma2=(parms[,5])^2,
			ylab="Change in EDD")
			
plot(outB,plotDif=TRUE)
}
}
\keyword{nonlinear}
\keyword{hplot}
