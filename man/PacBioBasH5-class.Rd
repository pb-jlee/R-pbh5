\name{PacBioBasH5-class}
\Rdversion{1.1}
\docType{class}
\alias{PacBioBasH5-class}
\alias{initialize,PacBioBasH5-method}
\alias{getQualityValue,PacBioBasH5-method}
\alias{getWidthInFrames,PacBioBasH5-method}
\alias{PacBioBasH5}

\title{Class "PacBioBasH5"}
\description{
  The base h5 (bas.h5) file represents basecalls. This file format is
  used to store information about called bases derived from pulses. This
  file-format additionally contains a subset of the pulse information
  used in the ascertainment of bases. It is a proper-subset of the
  information stored in pls.h5 files. 
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{PacBioBasH5(fileName)}.
}
\section{Slots}{
  \describe{
    \item{\code{baseEvents}:}{Object of class \code{"matrixOrNull"}}
    \item{\code{baseCallsG}:}{Object of class \code{"H5ObjOrNull"}}
    \item{\code{version}:}{Object of class \code{"stringOrNull"}}
    \item{\code{contents}:}{Object of class \code{"list"}}
    \item{\code{fileName}:}{Object of class \code{"character"}}
    \item{\code{ePtr}:}{Object of class \code{"externalptr"}}
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
    \item{initialize}{\code{signature(.Object = "PacBioBasH5")}:}
  }
}
\examples{
  showClass("PacBioBasH5")
}
\keyword{classes}
