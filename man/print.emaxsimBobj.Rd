\name{print.emaxsimBobj}
\alias{print.emaxsimBobj}
\title{
Print a summary of the fitted Emax model
}
\description{
Print a summary of the fitted Emax model.  
Printed output returned as invisible matrix.
}
\usage{
\method{print}{emaxsimBobj}(x, nprint=min(length(x$y),20), ...)
}
\arguments{
  \item{x}{
Object output by the extractor function [] for \code{\link{emaxsimB}}
}
\item{nprint}{ Number of observations to print.  If a vector of
length 2, \code{nprint} is the range of data to print. }
  \item{\dots}{
No options implemented.
}
}
