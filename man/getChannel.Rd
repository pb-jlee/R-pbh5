\name{getChannel}
\alias{getChannel}
\alias{getChi2}
\alias{getIsPulse}   
\alias{getMaxSignal}  
\alias{getMeanSignal}   
\alias{getMidSignal}  
\alias{getStartFrame}
\alias{getWidthInFrames}
\alias{getMidStdDevSignal}
\alias{getClassifierQV}
\title{
  Pulsecall-level Metrics 
}
\description{
  Retrieve matrices (or vectors) of features stored in the
  pls.h5 files. This information is not necessarily stored in the pls.h5
  file and code should check before depending on its existence. 
}
\usage{
  getChannel(plsH5, holeNumbers)
}
\arguments{
  \item{plsH5}{
    An object of class \code{PacBioPlsH5}
  }
  \item{holeNumbers}{
    Which ZMWs to grab.
  }
}
\value{
  A list of numeric vectors or matrices depending on the size of the
  dataset to be retrieved.
}
\details{
  In order to determine the meaning of the particular features consult
  the PacBio bas.h5 documentation.
}
\keyword{datasets}
