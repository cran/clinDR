\name{plotBdensity}
\alias{plotBdensity}

\title{Density plot displaying Bayes prior or posterior dose response}
\description{Density plot over a grid of doses displaying the prior 
or posterior distribution for the mean dose response computed from 
simulated input model parameters.}
  
\usage{
plotBdensity(dgrid, 
      parm, 
      modelFun=emaxfun,
      qlevL=c(0.025,0.05,0.10,0.25),
      plotDif= FALSE,
      logit= FALSE, ...)
}

\arguments{
  \item{dgrid}{The Bayes prior or posterior summaries are evaluated and plotted on
the \code{dgrid} dosing values }
  \item{parm}{Matrix of simultated parameter values (each row is a
simulated parameter vector).  The \code{parm} values
must be constructed for use in the model function \code{modFun}.  The
default is a 4-parameter Emax model with parameters
(log(ED50),lambda,Emax,E0).  For a 3-parameter model, set lambda=1 for each
simulated parameter vector.}
  \item{modelFun}{The mean model function.  The first argument is a
scalar dose, and the second argument is a matrix of parameter values.
The rows of the matrix are random draws of parameter vectors for the
model.  The default function is the 4-parameter Emax function \code{emaxfun}.}
  \item{qlevL}{Intervals are formed with percentile boundaries at (qlevL,1-qlevL).  
  \code{qlevL} must be increasing between (0,0.5).}
  \item{plotDif}{If TRUE, plot difference between doses and placebo.}
  \item{logit}{Default is F.  If T, inverse logit transform applied to
Emax function output for comparison to dose group sample proportions.}
  \item{\dots}{ Parameters passed to generic plot function }
}
\details{
  A sample of parameters from the joint prior or posterior distribution 
  must be supplied (typically produced by BUGS).  A density plot with
  contours corresponding to the perentiles in qlevL created by function
  \code{\link{DRDensityPlot}} . }

\value{
A list containing two matrices with the number of rows equal to the number dose grid points, and columns corresponding to percentiles in \code{qlevL}:
   \item{qL}{Lower perentiles from \code{qlevL} }
   \item{qH}{Upper percentiles 1-\code{qlevL}. }
}


\references{ Spiegelhalter, D., Thomas, A., Best, N., and Lunn, D. (2003),\emph{ WinBUGS User Manual Version 1.4}, Electronic version www.mrc-bsu.cam.ac.uk/bugs }
\author{ Neal Thomas }

\seealso{\code{\link{plot.plotB}}, \code{\link{plotD}},  \code{\link{plot.fitEmax}},
         \code{\link{DRDensityPlot}} }

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

dgrid<-seq(0,1,length=100)

plotBdensity(dgrid,parm=cbind(parms[,"led50"],parms[,"lambda"],
      parms[,"emax"],parms[,"e0[1]"]))

plotBdensity(dgrid,parm=cbind(parms[,"led50"],parms[,"lambda"],
      parms[,"emax"]),plotDif=TRUE,
       xlab='Dose',ylab='Dif with PBO')
}
\dontshow{
data("examples14")
exdat<-examples14[[6]]

parms<-matrix(c(-0.1665350, 0.3657811, -5.660137, 1.744753, 0.4050860,
-0.8463137, 0.3837361, -4.877676, 1.784098, 0.3943782,
-1.1811274, 0.3767222, -4.921861, 1.873861, 0.4266011,
 0.4729616, 0.3157714, -6.322768, 1.780517, 0.3646588,
 0.4255880, 0.3336959, -6.251558, 1.775438, 0.3657461),ncol=5,byrow=TRUE)

dgrid<-seq(0,1,length=5)

plotBdensity(dgrid,parm=parms[,1:4])
}
}
\keyword{nonlinear}
\keyword{hplot}
