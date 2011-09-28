\name{getChannel}

\alias{getChannel}
\alias{getChi2}
\alias{getIsPulse}   
\alias{getMaxSignal}  
\alias{getMeanSignal}   
\alias{getMidSignal}  
\alias{getStartFrame}
\alias{getWidthInFrames}

\title{
  Pulsecall-level Metrics 
}
\description{
  Retrieve matrices (or vectors) of kinetic features stored in the
  pls.h5 files. This information is not necessarily stored in the pls.h5
  file and code should check before depending on its existence. 
}
\usage{
  getChannel(plsH5, ...)
}
\arguments{
  \item{plsH5}{
    An object of class \code{PacBioPlsH5}
  }
  \item{\dots}{
    Extra arguments to the \code{getFromPlsH5} function. Notably,
    \code{holeNumbers} can be used to grab subsets of the ZMWs.
  }
}
\value{
  A list of numeric vectors or matrices depending on the size of the
  dataset to be retrieved.
}
\details{
  In order to determine the meaning of the particular features consult
  the XXX documentation.
}
\examples{
  plsH5 <- PacBioPlsH5(system.file("h5_files", "lambda_example.pls.h5", package = "pbh5"))
  mSignal <- Filter(Negate(is.null), getMaxSignal(plsH5))[[1]]
  head(mSignal)

  chi2 <- Filter(Negate(is.null), getChi2(plsH5))[[1]]
  boxplot(pchisq(chi2, df = 1))
}
\keyword{programming}
