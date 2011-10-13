\name{getPolymerizationRate}
\alias{getPolymerizationRate}
\title{
  Computes the Polymerization Rate for each Alignment
}
\description{
  'getPolymerizationRate' computes the polymerization rate for each
  alignment in the cmpH5 file.
}
\usage{
getPolymerizationRate(cmpH5, idx = seq.int(1, nrow(cmpH5)))
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to return.
  }
}
\details{
 Each alignment in the dataset has one polymerization rate value.
}
\value{
  Returns a vector of polymerization rates.
}

\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  pRate <- getPolymerizationRate(cmpH5)
  plot(density(pRate, na.rm = TRUE), xlab = "Polymerization Rate")
}
\keyword{datasets}

