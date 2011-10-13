\name{getCallsInRange}
\alias{getCallsInRange}
\alias{getConsensusInRange}
\alias{getCallsForReads}
\alias{getConsensusForReads}
\alias{computeConsensus}
\title{
  Function to get Calls In Range
}
\description{
  Returns a list containing a character vector for each position in the
  range refEnd - refStart + 1.
}
\usage{
getCallsInRange(cmpH5, refSeq, refStart = 1, refEnd =
                getRefLength(cmpH5, refSeq), idx = 1:nrow(cmpH5))

getConsensusInRange(cmpH5, refSeq, refStart = 1, refEnd = getRefLength(cmpH5, 
                    refSeq), idx = 1:nrow(cmpH5))

getCallsForReads(rawAlns, refStart, refEnd, strands, starts, ends) 

getConsensusForReads(rawAlns, refStart, refEnd, strands, starts, ends,
                     fast = TRUE) 
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{refSeq}{
  the reference sequence to retrieve.
}
  \item{refStart}{
  The reference start point.
}
  \item{refEnd}{
  The reference end point.
}
  \item{idx}{
    The indices of the alignments to retrieve.
  }
 \item{rawAlns}{
  The raw alignments to retrive with function 'getAlignmentsRaw.' 
}
 \item{strands}{
%%     ~~Describe \code{idx} here~~
}
 \item{starts}{
%%     ~~Describe \code{idx} here~~
}
 \item{ends}{
%%     ~~Describe \code{idx} here~~
}
 \item{fast}{
%%     ~~Describe \code{idx} here~~
}
}
\value{
  Functions 'getCallsInRange' and 'getCallsForReads'  return a lists of
  character vectors. Functions 'getConsensusInRange' and
  'getConsensusForReads' return character vectors.
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  callsInRange <- getCallsInRange(cmpH5, 1)
  head(callsInRange)
  cnsInRange <-getConsensusInRange(cmpH5, 1)
  head(cnsInRange)
  callsForReads <-getCallsForReads(getAlignmentsRaw(cmpH5, 1), 1, 2000, 1, 1, 2)
  head(callsForReads)
  cc <- computeConsensus(callsInRange)
  head(cc)
  cfr <- getConsensusForReads(getAlignmentsRaw(cmpH5, 1), 1, 2000, 1, 1, 2)
  head(cfr)
}
\keyword{datasets}

