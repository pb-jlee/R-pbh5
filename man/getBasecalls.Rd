\name{getBasecalls}
\alias{getBasecalls}
\alias{getDeletionTag}
\alias{getDeletionQV}
\alias{getInsertionQV}
\alias{getPulseIndex}
\alias{getPreBaseDeletionQV}
\alias{getSubstitutionTag}
\alias{getSubstitutionQV}
\alias{getPreBaseFrames}
\title{
  Get Basecall-level Metrics From pls.h5 or bas.h5 Files.
}
\description{
  These functions take an object of class bas.h5 or pls.h5 and return
  vectors for each ZMW represented by holeNumbers (all ZMWs) if no
  holeNumbers argument is given. Frequently, not all of these values
  are present in a bas.h5 or pls.h5 file and therefore users should
  either check or catch and error. 
}
\usage{
  getBasecalls(basH5, holeNumbers)
}
\arguments{
  \item{basH5}{
    An object of class \code{PacBioBasH5}.
  }
  \item{holeNumbers}{
    which ZMWs to grab.
  }
}
\value{
  A list of character vectors (or whatever object the function
  \code{convert} returns).
}
\examples{
basFiles <- list.files(system.file("h5_files", package = "pbh5"), pattern = "bas\\.h5", full.names = T)
cmpFiles <- list.files(system.file("h5_files", package = "pbh5"), pattern = "cmp\\.h5", full.names = T)

## instanteate the bas.h5 files and name them to work with some of the
## paired cmp.h5 bas.h5 functionality.
basH5s <- lapply(basFiles, PacBioBasH5)
names(basH5s) <- sapply(basH5s, getMovieName)

## returns the basecalls for holenumbers 1 and 2.
getBasecalls(basH5s[[1]], holeNumbers = 1:2)

## returns the most probable deleted base. 
getDeletionTag(basH5s[[1]], holeNumbers = 1:2)

## gets the index into the pulse calls.
getPulseIndex(basH5s[[1]], holeNumbers = 1:2)
}
\keyword{datasets}
