\name{readGroups}
\alias{readGroups}
\alias{refSeqs}
\title{
  Get readGroups and reference sequence names.
}
\description{
  Get Read group and Reference sequence names. These are used
  (esp. referenceSequence names) in different supporting functions such
  as \code{\link{getCoverageInRange}}
}
\usage{
readGroups(cmpH5)
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
}
\value{
  A character vector containing the names of the reference sequences or
  read groups.
}
\examples{
}
\keyword{programming}
