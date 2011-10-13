\name{applyZMWs}
\alias{applyZMWs}
\title{
  Apply Function to ZMWs
}
\description{
  Apply a function across ZMWs of a trace file. This function performs
  computation in blocks to avoid inevitable memory issues associated
  with the typical size of trc.h5 files.
}
\usage{
   applyZMWs(trcH5, fx, ..., .whichZMWs = NULL, .blockSize = 100)
}
\arguments{
  \item{trcH5}{An object of class trc.h5}
  \item{fx}{
    A function to apply to the channel by frame matrix.
  }
  \item{\dots}{
    Extra arguments to \code{fx}
  }
  \item{.whichZMWs}{
    To which ZMWs should the function be applied.
}
  \item{.blockSize}{
    Tuning parameter for adjusting the number of ZMWs that are brought
    into memory at any one time.
}
}
\details{
  This function is only mildly supported. The prospect of performing
  large scale computations on trc.h5 files in R is daunting.
}
\value{
  A list with the results of each application of \code{fx}
}
\seealso{
  \code{\link{PacBioTrcH5-class}}
}
\keyword{datasets}
