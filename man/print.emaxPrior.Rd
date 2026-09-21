\name{printemaxPrior}
\alias{print.emaxPrior}
\title{
  Generate and Save Protocol/SAP Text for Emax Prior Distributions
}
\description{
  Creates a templated description of the prior distribution for the
  parameters of the Emax model
  suitable for inclusion in a study protocol or statistical 
  analysis plan (SAP). The function outputs the text to a word file stored 
  by default in the current working directory.  The level of detail is 
  tailored to protocol or SAP requirements.
}
\usage{
  \method{print}{emaxPrior}(x, doc=FALSE, diffuse=c(pbo=TRUE,eff=TRUE), 
 file=paste('emaxPrior',format(Sys.time(), 
   "\%b\%d\%H\%M\%S"),'.docx',sep=''),
   modType=c('4','3'), docType=c('sap','protocol'),
   font.family = "Courier New", rounding=3, ...)
}
\arguments{
  \item{x}{Object created by function \code{emaxPrior.control}
  }
  \item{doc}{When \code{TRUE}, documentation for the prior distribution
    is returned.  Default is \code{FALSE} and the prior input is returned
    as a list for default printing}
  \item{diffuse}{When \code{TRUE}, the scale parameters are described
    as creating diffuse prior distributions for the corresponding 
    efficacy parameters.  When \code{FALSE}, sections are identified
    where the user must add justification for the informative 
    prior distributions.  
  }
  \item{file}{Filename to be stored in output location. If a custom filename
    is specified it must include the ".docx" extension.
    }
  \item{modType}{Character value ('4' or '3') that determines whether
    the 4-parameter sigmoidal Emax paramter is included, or the 3-parameter
    hyperbolic model is assumed with the Hill (slope) parameter set to 1
  }
  \item{docType}{When \code{'protocol'}, the prior description
    is less detailed.
  }
  \item{font.family}{The default 'Courier New' font family should be used 
   and changed only if it is not available.
  }
    \item{rounding}{Number of decimal places reported for numerical values.
  }
  \item{...}{No other inputs are supported.
  }
}
\details{
  If the emaxPrior object is entered at the command line, the implied print
  function is called without the required \code{diffuse} flag.  The
  object will be printed as a list without any generated text. 
  To generate a protocol/sap text, the \code{doc} must be \code{TRUE} and  
  the \code{diffuse} variable must be set (\code{TRUE/FALSE}).
}
\value{
  A Microsoft .docx file that can be edited for inclusion in
  a protocol/sap is written to the location specified by \code{file}. An
  invisible text version of the output is also returned.  If the output
  of \code{print.emaxPrior} is assigned, the text may then be printed 
  from the saved output using the function \code{cat}.
}
\author{
  Budhaditya Goswami,
  Neal Thomas
}
\seealso{
  \code{\link{emaxPrior.control}}
}
\examples{
  prior<-emaxPrior.control(epmu=0,epsca=4,difTargetmu=0,
                           difTargetsca=4,dTarget=20,
                           p50=(2+5)/2,
                           sigmalow=0.01,sigmaup=3)
  
  print(prior,doc=TRUE,diffuse=c(TRUE,TRUE),file=NULL)

}
\keyword{ prior }
\keyword{ Emax }
