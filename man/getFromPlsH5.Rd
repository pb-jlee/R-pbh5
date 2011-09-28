\name{getFromPlsH5}
\alias{getFromPlsH5}
\title{
  Get Datasets from pls.h5 Files
}
\description{
  Retrieve a dataset from a pls.h5 file which is stored as the base and
  pulse calls are stored. This function works with bas.h5 files as well,
  but only inasmuch as they contain the correct datasets.
}
\usage{
 getFromPlsH5(plsH5, plsH5GName = c("BaseCalls", "PulseCalls"), dsName, convert = NULL, holeNumbers = NULL)
}
\arguments{
  \item{plsH5}{
    An object of class \code{PacBioPlsH5} or \code{PacBioBasH5}
}
  \item{plsH5GName}{
    Either "BaseCalls" or "PulseCalls" depending on what group the
    dataset is contained within.
}
  \item{dsName}{
    The name of the dataset.
}
  \item{convert}{
    A function which can be applied to each slice of the dataset.
}
  \item{holeNumbers}{
    A vector of ZMW hole numbers representing the ZMWs one wishes to retrieve.
}
}
\details{
  Both pls.h5 and bas.h5 files store data associated with ZMWs in a
  structured format. This function implements grabbing according to that
  format. 
}
\value{
  A list of vectors or matrices.
}
\examples{
}
\keyword{programming}
