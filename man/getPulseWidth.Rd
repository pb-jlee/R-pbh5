\name{getPulseWidth}
\alias{getPulseWidth}
\title{
  Compute the Pulse Width
}
\description{
  'getPulseWidth' computes the pulsewidth for each base.
}
\usage{
getPulseWidth(cmpH5, idx)
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
  Pulsewidth is computed in seconds.
}
\value{
  'getPulseWidth' returns a list of vectors with pulsewidths.
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  pw <- getPulseWidth(cmpH5)
  logPW <- lapply(pw, function(x) log(x))
  plot(density(logPW[[1]], na.rm = TRUE), main = "Distribution of pulse
  widths from the first alignment")
}
\keyword{datasets}

