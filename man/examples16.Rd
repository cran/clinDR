\name{examples16}
\docType{data}
\alias{examples16}
\alias{meta16}
\title{Data from 2016 dose response meta-analysis}
\description{
Data from dose response meta-analysis published in 2016. 
}
\usage{data(examples16)
}
\format{
  The data for each of the 62 drugs in the analysis are stored in the list 
  \code{examples16}, which has length 62. Meta-data describing each drug 
  are stored in vectors of length 62 collected in a list \code{meta16}:  \code{drugs} (numerical id), 
  \code{drugname} (brand name),
  \code{ta} (therapuetic area), 
  \code{endpoint} (name of endpoint \code{y}),
  \code{epsource} (endpoint source, e.g., lab, PRO,...),
  \code{continuous} (type of primary outcome),
  \code{nprots}, (number of dose
  finding protocols), \code{ndoses} (number of doses studied).  The contents of 
  \code{examples16} are
  \describe{
    \item{\code{y}}{primary outcome, dose group mean or 0/1 for binary outcome}
    \item{\code{dose}}{total daily dose}
    \item{\code{prot}}{numerical protocol indicator}
    \item{\code{protnames}}{Sponsor-assigned protocol identifier}
    \item{\code{sem}}{SE for the y.  Computed using pooled SD for continuous
    data.  Based on binomial variance for binary data.}
    \item{\code{regimen}}{e.g., qd, bid, ...}
    \item{\code{nsize}}{sample size for each dose group/protocol}
     \item{\code{sd}}{within dose/protocol group pooled sd for continuous data}   
    \item{\code{esign}}{sign of intended effect}  
    \item{\code{pboAdj}}{\code{TRUE} when dose group means have placebo means 
    subtracted}  
    }
}

\details{
Details are given in the publication

Thomas, N., and Roy, D. (2016).  Analysis of clinical dose-response in small-molecule drug development:  2009-2014. Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317
<doi:10.1080/19466315.2016.1256229>}


\examples{
data(examples16)
objects(examples16[[1]])
names(meta16)
}

\keyword{datasets}
