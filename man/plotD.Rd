\name{plotD}
\alias{plotD}

\title{Basic plot of dose group means}
\description{
 Plot dose group means vs dose with options to connect points by
lines, and include CI about each dose group mean based on within-group
SDs
}
\usage{
plotD(y, dose, baseline, se = TRUE, line = TRUE, 
    meansOnly=FALSE,sem=NULL,clev = 0.8, 
    xlab='Dose',ylab='Response', logScale=FALSE)
}

\arguments{
  \item{y}{Outcomes }
  \item{dose}{Doses corresponding to outcomes }
  \item{baseline}{If present, ANACOVA means are plotted, adjusted for
baseline.  Baseline is optional.}
  \item{se}{If T, plot CI for each dose group.}
  \item{line}{If T, dose group means are connected by a line}
  \item{meansOnly}{If T, y contains dose group means rather than
  individual observations.  Baseline cannot be specified.}
  \item{sem}{If meansOnly and se=T, sem must contain the corresponding
  standard errors}
  \item{clev}{Level of CI for dose group means}
  \item{xlab}{Label for x-axis}
  \item{ylab}{Label for y-axis}
  \item{logScale}{If \code{TRUE}, log scale is used for dose.}
}

\value{
Returns a list with the ggplot object and two vectors with
the dose group means and their standard errors.} 

\author{Neal Thomas}

\seealso{ \code{\link{plot.fitEmax}}, \code{\link{plotB}} }
\examples{
data("metaData")
exdat<-metaData[metaData$taid==2 & metaData$etype==1,]
with(exdat,plotD(rslt,dose,meansOnly=TRUE,se=TRUE,sem=se,ylab=
"Y",xlab="Dose(mg)"))
}
\keyword{hplot}

