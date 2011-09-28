\name{getAccuracy}
\alias{getAccuracy}
\alias{getInsertionAndDeletionContext}
\alias{getReadLength}
\alias{getMachineName}
\alias{getPolymerizationRate}
\alias{getLocalPolymerizationRate}
\alias{getRefSeqNames}
\title{
  Compute Summary Statistics on Alignments
}
\description{
  These functions compute summary statistics for each alignment
  specified by \code{idx}. See Details for more information.
}
\usage{
getAccuracy(cmpH5, idx = seq.int(1, nrow(cmpH5)))

getLocalPolymerizationRate(cmpH5, idx = seq.int(1, nrow(cmpH5)), binFunction = function(aln) {
    cut(seq.int(1, nrow(aln)), 10)
})
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to return.
  }
  \item{binFunction}{
    A function which partitions each read into bins. Used in the local
    polymerization rate definition.
  }
}
\details{
  Accuracy is computed as: (nMismatches + nInsertions +
  nDeletions)/readLength.
}
\value{
  These functions typically return a single value for each alignment. An
    exception is the local polymerization rate calculation. This
    function returns a vector of values for each alignment. 
}
\examples{
  cmpH5 <- PacBioCmpH5(system.file("h5_files", "lambda_example_fromPls.cmp.h5", package = "pbh5"))
  rl <- getReadLength(cmpH5)
  
}
\keyword{programming}
