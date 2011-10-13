\name{getTemplateSpan}
\alias{getTemplateSpan}
\alias{getTemplateEnd}
\alias{getTemplateStart}
\alias{getTemplateStrand}
\alias{getReadLength}
\alias{getAccuracy}
\alias{getFullRefNames}
\alias{getRefNames}
\title{
  Functions to get information about the alignments.
}
\description{
  These functions retrieve information about alignments. Namely,
  position, length, accuracy.
}
\usage{
  getTemplateSpan(cmpH5, idx = seq.int(1, nrow(cmpH5)))
  getAccuracy(cmpH5, idx = seq.int(1, nrow(cmpH5)))
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to retrieve.
  }
}
\value{
  Return numeric vectors of length: \code{length(idx)}.
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5", package = "pbh5"))
  tStrand <- getTemplateStrand(cmpH5)
  tStart <- getTemplateStart(cmpH5)
  tEnd <- getTemplateEnd(cmpH5)
  tSpan <- getTemplateSpan(cmpH5)
  matrix <- cbind(tStart, tEnd, tSpan)
  head(matrix, 10)
}
\keyword{datasets}

