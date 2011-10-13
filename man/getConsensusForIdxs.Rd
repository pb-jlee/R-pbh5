\name{getConsensusForIdxs}
\alias{getConsensusForIdxs}
\title{
  Return Consensus For Indices
}
\description{
  Get the consensus sequence from the indices.
}
\usage{
  getConsensusForIdxs(cmpH5, idx = 1:nrow(cmpH5))
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to return.
  }
}
\value{
  Returns a vector of characters. 
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  head(getConsensusForIdxs(cmpH5))
}
\keyword{datasets}

