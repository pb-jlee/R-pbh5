\name{getConsensusCalls}
\alias{getConsensusCalls}
\alias{getConsensusConfidence}
\title{
  Get Consensus Calls
}
\description{
  Retrieve sequence consensus calls if available.
}
\usage{
getConsensusCalls(cmpH5, refSeq, collapse = FALSE)
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}
}
  \item{refSeq}{
    The reference sequence to retrieve.
  }
  \item{collapse}{
    Whether the consensus calls should be collapsed into a single
    character string or whether they should be returned as vector of
    characters. 
  }
}
\details{
  Consensus is obtained via the XXX program documented at: XXX
}
\value{
  Either a single character string or a vector of characters of length 1.
}
\references{
  XXX
}
\examples{
}
\keyword{programming}
