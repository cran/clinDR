\name{compileStanModels}
\alias{compileStanModels}
\title{ Compile \code{rstan} Emax models after package \code{clinDR} is installed

}
\description{
Compile \code{rstan} code for various Emax models used by \code{fitEmaxB}
and \code{emaxsimB}.  This function must be executed once after the 
\code{clinDR} package is installed.  
}
\usage{
compileStanModels()
}
\details{
The compiled models are stored in the \code{models} sub-directory of the installed
\code{clinDR} package.  The user must have write-access to the package directory.
The package can be installed in a user-specified directory if the user
does not have write privileges for the default 
package directory. Execution requires several minutes. The compiled models are
32- or 64- bit specific. Both sets must be compiled if the compiled R type
is changed (they are stored in sub-directories comp32 or comp64).  It is recommended
to execute the function again if the package \code{rstan} is updated.
}
\value{
Six model files with extensions '.rds' should be created in the package directory
in the sub-directory 'models'.
}
\author{
Neal Thomas
}

