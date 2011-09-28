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
  Alignments are encoded in a compressed format - more details of this
  format can be found at: 
}
\value{
  A list of integer vectors representing each alignment.
}
\examples{
}
\keyword{programming}