\name{getRefPath}
\alias{getRefPath}
\alias{getRefName}
\alias{getRefLength}
\title{
  Get data from refGroup table.
}
\description{
  Gets an element from refGroup table using any of name, index, or fasta
  name.
}
\usage{
  getRefPath(cmpH5, refSeq)
}
\arguments{
  \item{cmpH5}{
  }
  \item{refSeq}{Any of: index, name, or name from the fasta file.}
}
\value{
  character representing the absolute path to the reference group.
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5", package = "pbh5"))
  getRefPath(cmpH5, 1) == getRefPath(cmpH5, "ref000001")
}
\keyword{classes}

