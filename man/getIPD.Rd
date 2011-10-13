\name{getIPD}
\alias{getIPD}
\title{
  Compute the IPD
}
\description{
  'getIPD' computes the inter-pulse distance for each base.
}
\usage{
  getIPD(cmpH5, idx)
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
 Inter-pulse distances are given in seconds.
}
\value{
  A list of vectors of inter-pulse distances. 
}
\seealso{
  \code{\link{PacBioCmpH5-class}}
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  ipd <- getIPD(cmpH5)
  logIPD <- lapply(ipd, function(x) log(x))
  plot(density(logIPD[[1]], na.rm = TRUE), main = "Distribution of IPDs from the first alignment")
}
\keyword{datasets}
