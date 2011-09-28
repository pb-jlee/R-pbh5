\name{getClassifierQV}
\alias{getClassifierQV}
\title{
  Return Classifier Quality Values
}
\description{
  Return classifier quality values from a class/file which contains
  them. 
}
\usage{
  getClassifierQV(h5Obj, ...)
}
\arguments{
  \item{h5Obj}{
    Objects containing classifier QVs -- currently,
    \code{\link{PacBioCmpH5-class}} or \code{\link{PacBioPlsH5-class}} 
  }
  \item{\dots}{
    Extra arguments -- See class documentation for the particular \code{h5Obj}.
  }
}
\value{
  A list of vectors of classifier quality values.
}
\examples{
}
\keyword{programming}
