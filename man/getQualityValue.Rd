\name{getQualityValue}
\alias{getQualityValue}
\alias{getQualityValue,PacBioCmpH5-method}
\alias{getQualityValue,PacBioPlsH5-method}
\title{
  Return Classifier Quality Values
}
\description{
  Return classifier quality values from a class/file which contains
  them. 
}
\usage{
  getQualityValue(h5Obj, ...)
}
\arguments{
  \item{h5Obj}{
    An object of class \code{\link{PacBioCmpH5-class}} or
    \code{\link{PacBioBasH5-class}} 
  }
  \item{\dots}{
    Extra arguments -- See class documentation for \code{h5Obj}
  }
}
\value{
  A list of vectors of quality values.
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  values <- getQualityValue(cmpH5)
  class(values)
  length(values)
  length(values[[1]])
  head(values[[1]])
}
\keyword{datasets}
