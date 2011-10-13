\name{getBaseEvents}
\alias{getBaseEvents}
\alias{getPulseEvents}
\title{
  Get Base or Pulse Events
}
\description{
  Retrieve the base events matrix from either a bas.h5 or pls.h5
  file. This matrix contains the number of basecalls for each ZMW. This
  code is generally not used by the user, but rather is called by
  internal code for getting data out of [pls|bas].h5 files.
}
\usage{
  getBaseEvents(basH5)
}
\arguments{
  \item{basH5}{
   An object of class \code{PacBioBasH5}.
 }
}
\value{
  A matrix of base events.
}
\keyword{datasets}
