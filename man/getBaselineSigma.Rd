\name{getBaselineSigma}
\alias{getBaselineSigma}
\title{
  Get Baseline Sigma
}
\description{
  Retrieve the baseline sigma dataset from pls.h5 files.
}
\usage{
  getBaselineSigma(plsH5)
}
\arguments{
  \item{plsH5}{
    An object of class \code{PacBioPlsH5}
  }
}
\details{
  Baseline sigma is a 4-column matrix with nrow = number of ZMWs. Each
  value represents the channel-specific standard deviation.
}
\value{
  A 4-column matrix
}
\seealso{
  \code{\link{PacBioPlsH5-class}}
}
\keyword{datasets}

