\name{doWithPlsAndCmp}
\alias{doWithPlsAndCmp}
\alias{doWithPlsAndTrc}
\title{
  Do With Pulse And CpmH5 
}
\description{
  Executes function across a set of synced up PacBioCmpH5 and
  PacBioCmpH5 files.
}
\usage{
doWithPlsAndCmp(cmpH5, plsH5s, fx, ..., SIMPLIFY = TRUE)

doWithPlsAndTrc(plsH5, trcH5, holeNumbers, fx = rowMeans) 
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{plsH5s}{
    Objects of class \code{PacBioPlsH5}.
}
  \item{fx}{
    Statistic to compute, ie row means.
}
  \item{\dots}{
  Extra arguments
}
\item{SIMPLIFY}{
}
\item{plsH5}{
  An object of class \code{PacBioPlsH5}.    
}
\item{trcH5}{
  An object of class \code{PacBioTrcH5}.
}
\item{holeNumbers}{
}
}
\keyword{datasets}

