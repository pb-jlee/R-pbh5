\name{getCoverageInRange}
\alias{getCoverageInRange}
\alias{getReadsInRange}

\title{
  Get Coverage in Range
}
\description{
  Return a vector of length end - start + 1 of the number of reads
  overalapping that range. 
}
\usage{
getCoverageInRange(cmpH5, refSeq, start, end, idx = NULL)
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}
}
  \item{refSeq}{
    Which reference sequence to use.
}
  \item{start}{
    The start of the range (1-based inclusive).
  }
  \item{end}{
    The end of the range (1-based inclusive).
  }
  \item{idx}{
    The indices alignments from the alignmentIndex to use. This argument
    can be used to filter certain alignments.
  }
}
\value{
  A vector of length end - start + 1
}
\examples{
}
\keyword{programming}
