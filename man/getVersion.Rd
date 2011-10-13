\name{getVersion}
\alias{getVersion}
\title{
  Get Version Information from PacBio HDF5 files
}
\description{
  Retrieves any version information if it exists in the file.
}
\usage{
  getVersion(pacBioDataFile)
}
\arguments{
  \item{pacBioDataFile}{
    An object extending the class \code{PacBioDataFile}.
  }
}
\value{
  A character string or Null if no versioning information exists in the file.
}
\keyword{datasets}
