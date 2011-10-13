\name{PacBioAlnH5-class}
\Rdversion{1.1}
\docType{class}
\alias{PacBioAlnH5-class}
\alias{getAlignments,PacBioAlnH5-method}
\alias{nrow,PacBioAlnH5-method}
\alias{PacBioAlnH5}

\title{Class "PacBioAlnH5"}
\description{
  The align h5 (aln.h5) class represents pairwise alignments. This file
  is not produced by default and is one of the least well supported of
  the PacBio hdf5 file formats.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{PacBioAlnH5(fileName)}.
}
\section{Slots}{
  \describe{
    \item{\code{ZMWs}:}{Object of class \code{"matrix"} describing how
      the file can be indexed.}
    \item{\code{alignment}:}{Object of class \code{"H5Obj"} pointing to
      alignment group in hdf5 file.}
    \item{\code{version}:}{Object of class \code{"stringOrNull"}}
    \item{\code{fileName}:}{Object of class \code{"character"}}
    \item{\code{ePtr}:}{Object of class \code{"externalptr"} points to
      the H5File.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{PacBioDataFile}"}, directly.
  Class \code{"\linkS4class{H5File}"}, by class "PacBioDataFile", distance 2.
  Class \code{"\linkS4class{H5Obj}"}, by class "PacBioDataFile", distance 3.
  Class \code{"\linkS4class{H5ObjOrNull}"}, by class "PacBioDataFile", distance 4.
}
\section{Methods}{
  \describe{
    \item{getAlignments}{\code{signature(h5Obj = "PacBioAlnH5")}: ... }
    \item{nrow}{\code{signature(x = "PacBioAlnH5")}: ... }
  }
}
\examples{
  showClass("PacBioAlnH5")
}
\keyword{classes}
