\name{getReadsInRange}
\alias{getReadsInRange}

\title{
  Retrieves Reads in Specified Range
}
\description{
  'getReadsInRange' retrieves reads in range based on arguments refSeq,
  refStart, refEnd and idx.
}
\usage{
getReadsInRange(cmpH5, refSeq, refStart = 1, refEnd = NA, idx = 1:nrow(cmpH5))
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{refSeq}{
  Reference sequence from which to retrieve reads.
}
  \item{refStart}{
  Reference start point from which to retrieve reads.
}
  \item{refEnd}{
  Reference end point from which to retrieve reads.
}
  \item{idx}{
    The indices of the alignments to return.
  }
}
\value{
  'getReadsInRange' returns a vector of integers.
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  readsInRange <- getReadsInRange(cmpH5, 1, 100, 200)
  class(readsInRange)
  head(readsInRange)

}
\keyword{datasets}

