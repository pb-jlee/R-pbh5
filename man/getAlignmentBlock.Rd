\name{getAlignmentBlock}
\alias{getAlignmentBlock}
\title{
  Returns a Specified Alignment Block
}
\description{
  Returns a view of the MSA for a group of reads.
}
\usage{
  getAlignmentBlock(cmpH5, ref, refStart, refEnd)
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{ref}{
    Reference sequence
  }
  \item{refStart}{
    Reference start point.
  }
  \item{refEnd}{
    Reference end point.
  }
}
\value{
  Returns a matrix of characters.
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  getAlignmentBlock(cmpH5, 1, 20, 40)
}
\keyword{datasets}

