\name{getCoverageInRange}
\alias{getCoverageInRange}
\title{
  Get Coverage in Range
}
\description{
  Return a vector of length end - start + 1 of the number of reads
  overalapping that range. 
}
\usage{
 getCoverageInRange(cmpH5, refSeq, refStart, refEnd, idx = NULL)
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}
}
  \item{refSeq}{
    Which reference sequence to use.
}
  \item{refStart}{
    The start of the range (1-based inclusive).
  }
  \item{refEnd}{
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
  require(pbh5)  
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  coverage <- getCoverageInRange(cmpH5, 1, 1000, 2000)
  boxplot(coverage, ylab = "Pileup Coverage") 
}
\keyword{datasets}
