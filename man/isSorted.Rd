\name{isSorted}
\alias{isSorted}
\title{
  Find Out Whether cmpH5 Object is Sorted
}
\description{
  'isSorted' returns a boolean statement whether the cmpH5 object is sorted.
}
\usage{
isSorted(cmpH5)
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
}
\value{
  Returns a logical value.
}
%% ~~ Need working examples here. ~~
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5")) 
}
\keyword{datasets}

