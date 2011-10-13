\name{PacBioPlsH5-class}
\Rdversion{1.1}
\docType{class}
\alias{PacBioPlsH5-class}
\alias{getSNR,PacBioPlsH5-method}
\alias{getClassifierQV,PacBioPlsH5-method}
\alias{getWidthInFrames,PacBioPlsH5-method}
\alias{PacBioPlsH5}

\title{Class "PacBioPlsH5"}
\description{
  The pulse h5 (pls.h5) class represents basecalls and
  pulsecalls. Pulsecalls are translated to basecalls and that mapping is
  contained in the pls.h5 file. The pulse h5 class "contains" an bas.h5
  file. 
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{PacBioPlsH5(fileName)}.
}
\section{Slots}{
  \describe{
    \item{\code{pulseEvents}:}{Object of class \code{"matrixOrNull"}}
    \item{\code{pulseCallsG}:}{Object of class \code{"H5ObjOrNull"}}
    \item{\code{baseEvents}:}{Object of class \code{"matrixOrNull"}}
    \item{\code{baseCallsG}:}{Object of class \code{"H5ObjOrNull"}}
    \item{\code{version}:}{Object of class \code{"stringOrNull"}}
    \item{\code{fileName}:}{Object of class \code{"character"}}
    \item{\code{ePtr}:}{Object of class \code{"externalptr"}}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{PacBioBasH5}"}, directly.
  Class \code{"\linkS4class{PacBioDataFile}"}, by class "PacBioBasH5", distance 2.
  Class \code{"\linkS4class{H5File}"}, by class "PacBioBasH5", distance 3.
  Class \code{"\linkS4class{H5Obj}"}, by class "PacBioBasH5", distance 4.
  Class \code{"\linkS4class{H5ObjOrNull}"}, by class "PacBioBasH5", distance 5.
}
\section{Methods}{
  \describe{
    \item{getQualityValue}{\code{signature(h5Obj = "PacBioPlsH5")}: ... }
  }
}
\examples{
showClass("PacBioPlsH5")
}
\keyword{classes}
