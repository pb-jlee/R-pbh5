\name{getIPD}
\alias{getIPD}
\alias{getLight}
\alias{getPulseWidth}
\alias{getStartTime}
\alias{getPkmid}

\title{
  Get Pulse and Basecall Features from cmp.h5 Files
}
\description{
  Retrieve datasets stored within cmp.h5 files. These datasets might not be available
  in all cmp.h5 files. 
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
  Pulse and basecall information are often stored in cmp.h5 files. This
  information can be accessed using accessor functions. The meaning of
  this information can be found in the XXX documentation.
}
\value{
  A list of vectors of the particular dataset which was retrieved. 
}
\seealso{
  \code{\link{PacBioCmpH5-class}}
}
\examples{
   cmpH5 <- PacBioCmpH5(system.file("h5_files", "lambda_example_fromPls.cmp.h5", package = "pbh5"))
   plot(density(getLight(cmpH5)[[1]], na.rm = TRUE))
}
\keyword{programming}
