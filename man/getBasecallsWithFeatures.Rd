\name{getBasecallsWithFeatures}
\alias{getBasecallsWithFeatures}
\title{
  Get Bascalls with Pulse Features
}
\description{
  Retrieve basecalls with attendant pulse features. These features must
  be present in the pls.h5 file.
}
\usage{
getBasecallsWithFeatures(plsH5, features = .ofPulses, holeNumbers, collapse = FALSE)
}
\arguments{
  \item{plsH5}{
    An object of class \code{PacBioPlsH5}.
  }
  \item{features}{
    A character vector of feature names. 
  }
  \item{holeNumbers}{
    A vector of ZMW numbers. 
  }
  \item{collapse}{
    Whether a list of data.frames or a single data.frame should be returned.
  }
}
\details{
  This function returns pulse features and basecalls. The available
  pulse features are described in \code{\link{PacBioPlsH5-class}}.
}
\value{
  A list of data.frames with basecalls and pulse features.
}
\seealso{
  \code{\link{PacBioCmpH5-class}}, \code{\link{PacBioBasH5-class}}
}
\examples{
}
\keyword{programming}
