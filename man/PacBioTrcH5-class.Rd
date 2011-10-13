\name{PacBioTrcH5-class}
\Rdversion{1.1}
\docType{class}
\alias{PacBioTrcH5-class}
\alias{PacBioTrcH5}
\alias{makePlsAndTrcProcessor}

\title{Class "PacBioTrcH5"}
\description{
  The trace h5 (trc.h5) class represents trace information. Trace files
  are *massive* and therefore they are often deleted immediately. 
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{PacBioTrcH5(fileName, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{version}:}{Object of class \code{"stringOrNull"} ~~ }
    \item{\code{fileName}:}{Object of class \code{"character"} ~~ }
    \item{\code{ePtr}:}{Object of class \code{"externalptr"} ~~ }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{PacBioDataFile}"}, directly.
  Class \code{"\linkS4class{H5File}"}, by class "PacBioDataFile", distance 2.
  Class \code{"\linkS4class{H5Obj}"}, by class "PacBioDataFile", distance 3.
  Class \code{"\linkS4class{H5ObjOrNull}"}, by class "PacBioDataFile", distance 4.
}
\section{Methods}{
  No methods defined with class "PacBioTrcH5" in the signature.
}
\examples{
showClass("PacBioTrcH5")
}
\keyword{classes}
