\name{getAlignmentsWithFeatures}
\alias{getAlignmentsWithFeatures}
\title{
  Get Alignments with Features
}
\description{
  Retrieve alignments with attendant pulse features.
}
\usage{
getAlignmentsWithFeatures(cmpH5, features, idx = seq.int(1, nrow(cmpH5)), collapse = FALSE)
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{features}{
    A character vector naming the additional features.
  }
  \item{idx}{
    The indices of the alignments to return.
  }
  \item{collapse}{
    Whether a list of data.frames should be returned or a single
    data.frame with a column naming the alignment.
  }
}
\value{
  Either a list of data.frames or a data.frame
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  alnsPlusFeatures <- getAlignmentsWithFeatures(cmpH5,
        features = c("IPD","PulseWidth", "pkmid"), idx = 1, collapse = TRUE)
  head(alnsPlusFeatures)
  boxplot(IPD ~ reference, data = alnsPlusFeatures, log = 'y',
           ylim = c(.05, 10), main = "Distribution of IPDs Conditioned on Base")
  plot(density(log(alnsPlusFeatures$PulseWidth), na.rm = TRUE), main =
 "Distribution of Global Pulsewidth")
}
\keyword{datasets}
