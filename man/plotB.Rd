\name{plotB}
\alias{plotB}

\title{Plot Bayes  dose response curve and dose group means  }
\description{Plot a dose response curve fit by Bayes MCMC methods
(with optional posterior interval bars).  Also plot dose
group means (with optional CI bars)

}
\usage{
plotB(y, 
	dose, 
	parm, 
	sigma2,
	count=rep(1,length(y)),
	dgrid=sort(unique(c(seq(0,max(dose),length=50), dose))),
	predict= TRUE,plotDif=FALSE,plotMed=FALSE,
	plotResid=FALSE,clev=0.8,
	binary=c('no','logit','probit','BinRes'),BinResLev,
	BinResDir=c('>','<'),
	activeControl=FALSE,ac,yac,
	countac=rep(1,length(yac)),
	labac='Act Comp',shapeac=8,colac='red',
	symbol,symbolLabel='Group',symbolShape=8,
	symbolColor='red',symbolSize=4,
	xlim,ylim,xat=NULL,xlab="Dose",
	ylab=ifelse(plotDif,"Diff with Comparator","Mean"),
	modelFun=emaxfun,makePlot=TRUE,
	...)
}

\arguments{
  \item{y}{Outcomes, which may be sample means (see \code{counts}).  LSmeans from a saturated anacova 
  model can be supplied, in 
  which case it is assumed that the Bayesian dose response model also included the additive
  baseline covariates.}
  \item{dose}{Doses corresponding to outcomes }
  \item{parm}{Matrix of simultated parameter values (each row is a
simulated parameter vector).  The \code{parm} values
must be constructed for use in the model function \code{modFun}.  The
default is a 4-parameter Emax model with parameters
(log(ED50),lambda,Emax,E0).  For a 3-parameter model, set lambda=1 for each
simulated parameter vector.}
  \item{sigma2}{Simulated draws from the residual variance (assumed
additive, homogeneous).  The length of \code{sigma2} must be the same
as the number of rows of \code{parm}.  Set \code{sigma2} to all ones
for binary data.}
  \item{count}{Sample sizes for means-only summarized data.}
  \item{dgrid}{The Bayes posterior summaries are evaluated and plotted on
the \code{dgrid} dosing values }
  \item{predict}{If TRUE(default), the plotted intervals are
predictive intervals for the dose group sample means.}
  \item{plotDif}{Plot difference between doses and placebo. 
  It is assumed the lowest dose is placebo. If \code{activeControl}, the
  difference is with the active control mean, and the active 
  controls are not plotted.}
  \item{plotMed}{If TRUE, model-based curves are medians rather than means.}
  \item{plotResid}{If TRUE, a plot of the residuals formed from the
  dose group means minus the posterior dose group means.}
  \item{clev}{Level for confidence and Bayes intervals }
  \item{binary}{If \code{binary} is 'logit' or 'probit', \code{y} 
  is assumed to be binary and the appropriate backtransformation is applied
  to the Emax model output.  If \code{binary} is 'BinRes', the continuous
  variable \code{y} is converted to a binary responder variable
  using \code{BinResLev} and \code{BinResDir}.  The 
  continuous Emax model output is converted to binary estimation and
  prediction assuming normally distributed residuals.}
   \item{BinResLev}{A cut level for a responder variable formed from a 
   continuous endpoint.  Rates are computed from the (continuous
   outcome) model parameters assuming normally distributed residuals.  The
   input \code{y} variable is converted to a responder variable.}
   \item{BinResDir}{If BinResDir='>', the responder variable is 1 when
   \code{y} is greater than the cut level, otherwise, it is 1 when \code{y} is less
   than the cut level.}
  \item{activeControl}{When TRUE, active comparator data must be supplied.  Each
  dose group (including PBO) are compared to the active comparator rather than PBO.}
  \item{ac}{Simulations from the posterior distribution of the mean response on 
  active comparator.  The number of simulations must match those for the dose
  response model.  For binary data, the simulated values must be transformed to the
  proportion scale.  This differs from the simulated model parameters.}
  \item{yac}{Outcomes for the active comparator group.  The coding conventions for 
  \code{y} are used.}
  \item{countac}{Sample sizes for summarized data corresponding to \code{count}.}
  \item{labac}{x-axis label for the active control group.}
  \item{shapeac}{Shape of the symbol for the active control group.}
  \item{colac}{Color of the symbol for the active control group.}
  \item{symbol}{An optional grouping variable for the dose group sample means.}  
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
  \item{xlim}{Plot limits for the x-axis}
  \item{ylim}{Plot limits for the y-axis}
    \item{xat}{The points at which tick-marks are to be drawn. Errors occur if the 
  points are outside the range of xlim. By default (when NULL) tickmark locations are    computed.}
  \item{xlab}{x-axis label }
  \item{ylab}{y-axis label }
  \item{modelFun}{The mean model function.  The first argument is a
scalar dose, and the second argument is a matrix of parameter values.
The rows of the matrix are random draws of parameter vectors for the
model.  The default function is the 4-parameter Emax function \code{emaxfun}.}
  \item{makePlot}{If FALSE, create numerical output but no plot.}
  \item{\dots}{ Parameters passed to generic plot function (not used) }
}
\details{
  A sample of parameters from the joint posterior distribution must be supplied (typically produced by BUGS).
The Bayesian dose response curve is the Bayes posterior mean (or median) at each
value on \code{dgrid}.  The bar (interval) is the (clev/2,1-clev/2)
Bayes posterior interval (which can differ from the Bayes HPD
interval).  The intervals are plotted only at the dose levels included
in the study.  Predictive intervals  are formed by
adding independent random draws from the sampling distributions of the dose group
sample means to the population means.

The function generates random numbers when \code{predict=TRUE}, so the random number generator/seed must
be set before the function is called for exact reproducibility.
}
\value{
Returns an object of class \code{plotB}.  Three  inputs are saved for
later plotting:  doses in the original design, dgrid, and clev.  
The following matrices are saved:

\item{pairwise}{The dose group means and their 
differences with placebo.  If a \code{baseline} is
supplied, the means are lsmeans adjusted to the mean baseline value.}
\item{modelABS}{Model-based posterior mean, median, 
posterior (clev/2,1-clev/2) intervals for the population means and
sample means.  One row per dose group}
\item{modelABSG}{Same as \code{modelABS} but computed on the input
grid of doses.}
\item{modelDIF}{Same as modelABS but with differences from placebo.}
\item{modelDIFG}{Same as modelDIF  but computed on the input
grid of doses.}

}

\note{
PlotB can also be used with draws from a prior distribution to
evaluate the prior dose response curve.
}

\references{ Spiegelhalter, D., Thomas, A., Best, N., and Lunn, D. (2003),\emph{ WinBUGS User Manual Version 1.4}, Electronic version www.mrc-bsu.cam.ac.uk/bugs }
\author{ Neal Thomas }

\seealso{\code{\link{plot.plotB}}, \code{\link{plotD}},  \code{\link{plot.fitEmax}} }
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
