\name{mcmc.control}
\alias{mcmc.control}
\title{
Settings for \code{restan} execution in function \code{fitEmaxB}
}
\description{
Set MCMC controls.  Also control spread of initial parameter values.
}
\usage{
mcmc.control(chains = 1, thin = 1, 
warmup = 500, iter = 5000* thin, 
propInit = 0.25, seed = 12357, adapt_delta = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chains}{
Number of chains
}
  \item{thin}{
Number of discarded sampled parameter values.  \code{warmup} and \code{iter} include
\code{thin}, so for example, to output \code{1000} samples, \code{iter} must be 
\code{1000} times \code{thin}.
}
  \item{warmup}{
See \code{rstan} documentation for function \code{sampling}.   
}
  \item{iter}{
See \code{rstan} documentation for function \code{sampling}.  
}
  \item{propInit}{
Initial values for \code{E0} and \code{Emax} are derived from the prior mean plus/minus
\code{propInit} times the prior SD.  \code{propInit} can be set to a small proportion if
very diffuse prior distributions are specified.
}
  \item{seed}{
Seed passed to \code{rstan}.
}
  \item{adapt_delta}{
See \code{rstan} documentation for function \code{sampling}.
}
}

