\name{getAlignments}
\alias{getAlignments}
\alias{templatePositionFromAln}
\alias{getAlignments,PacBioCmpH5-method}
\title{
  Return Alignments
}
\description{
  Return alignments from a class which contains alignments.
}
\usage{
  getAlignments(h5Obj, ...)
}
\arguments{
  \item{h5Obj}{
    An object containing alignments.
  }
  \item{\dots}{
    Extra arguments -- See the particular class documentation for
    specific methods.
  }
}
\details{
  Currently, alignments are two column matrices of mode character. 
}
\value{
  A list of alignment matrices.
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  nucleotideSeq <- getAlignments(cmpH5)
  class(nucleotideSeq)
  head(nucleotideSeq[[1]])
}
\keyword{datasets}
