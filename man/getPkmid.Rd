\name{getPkmid}
\alias{getPkmid}
\title{
  Computes the Pkmid Values
}
\description{
  'getPkmid' computes the Pkmid values, a component of Nignal to Noise
  Ratio (SNR). A vector of Pkmid values is computed for each alignment
  in the cmph5 file. 
}
\usage{
getPkmid(cmpH5, idx)
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
  Pkmid values are calculated from Pulse Sigma and Baseline Sigma values.
}
\value{
  'getPkmid' returns a list of vectors with numeric values. 
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  pkmid <- getPkmid(cmpH5, idx = 1:10)
  boxplot(pkmid,
  main = "Distribution of Pkmid values for the first ten alignments")
}
\keyword{datasets}

