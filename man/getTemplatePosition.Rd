\name{getTemplatePosition}
\alias{getTemplatePosition}
\title{
  Get Position in Template for Aligned Reads
}
\description{
  Returns a vector of template positions for each alignment
  positions. Essentially, it provides a mapping between alignment
  position to template position.
}
\usage{
  getTemplatePosition(cmpH5, idx, withAlignments = TRUE)
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to retrieve.
}
  \item{withAlignments}{
    Should the alignment also be returned.
  }
}
\value{
  A list of numeric vectors or a list of data.frames with the template
  position and alignment stored together.
}

\examples{
}
\keyword{programming}
