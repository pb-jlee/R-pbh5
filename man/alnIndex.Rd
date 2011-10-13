\name{alnIndex}
\alias{alnIndex}
\title{
  Get the AlignmentIndex from a PacBioCmpH5 object
}
\description{
  'alnIndex' returns the alignment index from a PacBioCmpH5 object, including
  information about ID, alnGroupPath, movieName, refName, tStart, tEnd,
  alignedStrand, holeNumber, setNumber, strobeNumber, moleculeID,
  rStart, rEnd, mapQV, nMatches, nMisMatches, nInsertions, ndeletions,
  offsetBegin, offsetEnd, nBackRead and  nOverlap.
}
\usage{
  alnIndex(cmpH5)
}
\arguments{
  \item{cmpH5}{An object of class PacBioCmpH5}
}
\value{
  'alnIndex' returns a data.frame where each row represents an
  alignment, and each column represents unique information about the
  alignment. 'alnGroup' returns a data.frame where each row is a unique
  Path. See description.  
}
\seealso{
  \code{\link{PacBioCmpH5-class}}
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  alignments <- alnIndex(cmpH5)
  class(alignments)
  head(alignments)
  alignGroup <-alnGroup(cmpH5)
  head(alignGroup)
}
\keyword{datasets}
