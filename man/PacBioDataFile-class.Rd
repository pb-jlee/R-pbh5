\name{PacBioDataFile-class}
\Rdversion{1.1}
\docType{class}
\alias{PacBioDataFile-class}
\alias{initialize,PacBioDataFile-method}
\alias{show,PacBioDataFile-method}

\title{Class "PacBioDataFile"}
\description{
  The parent class of all pbh5 classes. 
}
\section{Objects from the Class}{
  Objects of this class should not be directly instantiated.
}
\section{Slots}{
  \describe{
    \item{\code{version}:}{Object of class \code{"stringOrNull"} ~~ }
    \item{\code{fileName}:}{Object of class \code{"character"} ~~ }
    \item{\code{ePtr}:}{Object of class \code{"externalptr"} ~~ }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{H5File}"}, directly.
  Class \code{"\linkS4class{H5Obj}"}, by class "H5File", distance 2.
  Class \code{"\linkS4class{H5ObjOrNull}"}, by class "H5File", distance 3.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "PacBioDataFile")}: ... }
    \item{show}{\code{signature(object = "PacBioDataFile")}: ... }
  }
}
\examples{
showClass("PacBioDataFile")
}
\keyword{classes}
