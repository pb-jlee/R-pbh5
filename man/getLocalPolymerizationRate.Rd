\name{getLocalPolymerizationRate}
\alias{getLocalPolymerizationRate}
\title{
  Compute the local Polymerization Rate
}
\description{
  'getLocalPolymerizationRate' computes the local polymerization rate
  within each alignment by dividing each alignment into 10 sections.
}
\usage{
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
    A function which partitions each alignment into bins. The default
    option is 10 bins per alignment.
}
}
\value{
  'getLocalPolymerizationRate' returns a list of arrays with start and
  end times as well as corresponding local polymerization rates within
  those time-blocks.

}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  localRate <- getLocalPolymerizationRate(cmpH5, idx = 1:10)
  boxplot(localRate, ylab = "Local Polymerization Rate", xlab = "Index Numbers")
}
\keyword{datasets}

