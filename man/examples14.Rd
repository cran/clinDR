\name{examples14}
\docType{data}
\alias{examples14}
\alias{meta14}
\title{Data from 2014 dose response meta-analysis}
\description{
Data from dose response meta-analysis published in 2014. 
}
\usage{data(examples14)
}
\format{
  The data for each of the 37 drugs in the analysis are stored in the list 
  \code{examples14}, which has length 37. Meta-data describing each drug 
  are stored in vectors of length 37 collected in a list \code{meta14}:  \code{drugs} (numerical id), 
  \code{ta} (therapuetic area), 
  \code{endpoint} (name of endpoint \code{y}),
  \code{epsource} (endpoint source, e.g., lab, PRO,...),
  \code{continuous} (type of primary outcome),
  \code{nprots}, (number of dose
  finding protocols), \code{ndoses} (number of doses studied).  The contents of 
  \code{examples14} are
  \describe{
    \item{\code{y}}{primary outcome, dose group mean or 0/1 for binary outcome}
    \item{\code{dose}}{total daily dose}
    \item{\code{prot}}{numerical protocol indicator}
    \item{\code{sem}}{SE for the y.  Computed using pooled SD for continuous
    data.  Based on binomial variance for binary data.}
    \item{\code{regimen}}{e.g., qd, bid, ...}
    \item{\code{nsize}}{sample size for each dose group/protocol}
     \item{\code{sd}}{within dose/protocol group pooled sd for continuous data}   
    }
}

\details{
Details are given in the publication

Thomas, N., Sweeney, K., and Somayaji, V. (2014).  Meta-analysis of clinical dose response in a large drug development portfolio, Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317.
<doi:10.1080/19466315.2014.924876>

}

\examples{
data(examples14)
objects(examples14[[1]])
names(meta14)
}

\keyword{datasets}
