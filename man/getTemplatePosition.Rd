\name{getTemplatePosition}
\alias{getTemplatePosition}
\alias{getByTemplatePosition}
\title{
  Get Position in Template for Aligned Reads
}
\description{
  'getTemplatePosition' returns a vector of template positions for each
  alignment positions. Essentially, it provides a mapping between alignment
  position to template position. 'getBytemplatePosition' returns a
  data.frame with additional information on a position basis (default
  argument is 'getIPD.'
}
\usage{
  getTemplatePosition(cmpH5, idx, withAlignments = FALSE, asDataFrame = FALSE)
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
  \item{asDataFrame}{
    Whether to turn the results into a data.frame. 
  }  
}
\value{
  A list of numeric vectors or a list of data.frames with the template
  position and alignment stored together.
}

\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  tpos <- getTemplatePosition(cmpH5, 1, TRUE)
  tposIPD <- getByTemplatePosition(cmpH5, idx = 1)
  head(tposIPD)
  boxplot(elt ~ ref, ylim = c(0,5), data = tposIPD)
}
\keyword{datasets}
