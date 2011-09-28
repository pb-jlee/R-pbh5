\name{getBasecall}

\alias{getBasecall}
\alias{getDeletionTag}
\alias{getDeletionQV}
\alias{getInsertionQV}
\alias{getPulseIndex}
\alias{getPreBaseDeletionQV}
\alias{getSubstitutionTag}
\alias{getSubstitutionQV} 

\title{
  Get Basecall-level Metrics From pls.h5 Files.
}
\description{
  These functions take an object of class bas.h5 and return vectors for
  each ZMW represented by holeNumbers (all ZMWs) if no holeNumbers
  argument is given. Frequently, not all of these values are present in
  a bas.h5 or pls.h5 file and therefore users should either check or
  catch and error. 
}
\usage{
  getBasecall(basH5, convert = .f, ...)
}
\arguments{
  \item{basH5}{
    An object of class \code{PacBioBasH5}.
  }
  \item{convert}{
    The function to convert the underlying vectors returned. Each vector
  represents the results (pulsecalls or basecalls) from a ZMW. The
  canonical example is to convert ascii digits to characters. 
  }
  \item{\dots}{
    Extra arguments, of which, \code{holeNumbers} is most relevant. This
    argument can be used to specify which ZMWs to grab.
  }
}
\value{
  A list of character vectors (or whatever object the function
  \code{convert} returns).
}
\seealso{
  \code{\link{getFromPlsH5}}
}
\examples{
basH5 <- PacBioBasH5(system.file("h5_files", "lambda_example.bas.h5", package = "pbh5"))
basecalls <- getBasecall(basH5)
summary(sapply(basecalls, length))
}
\keyword{programming}
