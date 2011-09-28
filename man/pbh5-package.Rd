\name{pbh5-package}
\alias{pbh5-package}
\alias{pbh5}
\docType{package}
\title{
  Interact with Pacific Biosciences HDF5 Files
}
\description{
  The pbh5 package provides and interface to the various HDF5 files
  produced by Pacific Biosciences software. These files can be used to
  access unique features of PacBio data not readily available in
  standard formats.

  The package is structured around a set of hdf5 file formats and the
  classes used to represent them. The hdf5 format is a robust
  infrastructure for the construction of binary file formats. More
  information can be found at: www.hdfgroup.org. Roughly, it is a
  filesystem for your data. At the R level, the pbh5 package depends on
  the h5r package which is an open-source package located on
  r-forge. This package provides a generic read-only interface to hdf5
  files.

  The core classes defined in the pbh5 package are: PacBioCmpH5
  (pronounced comp h5), PacBioBasH5 (pronounced base h5), PacBioPlsH5
  (pronounced pulse h5). The pbh5 package also supports two other
  formats which will be encountered less frequently, if at all:
  PacBioTrcH5 (pronounced trace h5) and PacBioAlnH5 (pronounced align
  h5).

  All classes should be considered read-only in the sense that any
  modification of the underlying properties will have undefined
  consequences. Additionally, access to the data stored in the classes
  should proceed (almost) entirely through functions and methods defined
  herein. To the extent possible, objects are returned as simple R
  objects.

  Finally, and most importantly, both the code and the documentation are
  works-in-progress. To the extent that example files and code disagrees
  with the specification, the specifications should be trusted and a bug
  should be reported to <jbullard@pacificbiosciences.com>.
}
\details{
  \tabular{ll}{
    Package: \tab pbh5\cr
    Version: \tab 1.0\cr
    Date: \tab 2010-07-02\cr
    License: \tab BSD \cr
    LazyLoad: \tab yes\cr
    Depends: \tab methods, h5r\cr
  }
}
\author{
  Maintainer: <jbullard@pacificbiosciences.com>
}
\keyword{ package }
\seealso{
 The h5r package. 
}
\examples{
 ls('package:pbh5')
}
