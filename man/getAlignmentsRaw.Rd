\name{getAlignmentsRaw}
\alias{getAlignmentsRaw}
\title{
  Return Encoded Alignments
}
\description{
  Returns the encoded alignments as described in the specification for
  cmp.h5 files.
}
\usage{
getAlignmentsRaw(cmpH5, idx)
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to retrieve.
  }
}
\details{
  Alignments are encoded in a compressed format.
}
\value{
  A list of integer vectors representing each alignment.
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  rawAlns <- getAlignmentsRaw(cmpH5)
  head(rawAlns[[1]])
}
\keyword{datasets}