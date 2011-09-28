\name{associateCmpWithPls}
\alias{associateCmpWithPls}
\title{
  Associate cmp.h5 data with pls.h5 data
}
\description{
  Associate pulse data from pls.h5 files to alignments in a cmp.h5
  file. This function is decidedly slow and should be avoided for
  anything but exploratory analysis.
}
\usage{
associateCmpWithPls(cmpH5, plsFiles, features, idx = seq.int(1, nrow(cmpH5)), collapse = TRUE)
}
\arguments{
  \item{cmpH5}{
    An object of class PacBioCmpH5
  }
  \item{plsFiles}{
    A list of objects of class PacBioPlsH5 -- the movies represented by
    these files must be those stored in the \code{cmpH5} object.
}
  \item{features}{
    The features which should be associated.
}
  \item{idx}{
    The indices of the alignmentIndex to be joined.
}
  \item{collapse}{
    Whether the resulting object should be a list or a collapsed data.frame.
  }
}

\value{
  Either a data.frame or list of data.frames with the associated
  alignments and pulse features.
}
\details{
  cmp.h5 files can optionally contain pulse information derived from
  their parent pls.h5 file(s). This information is best copied to the
  cmp.h5 file in batch during the day-to-day production of these
  files. If, however, this does not happen by default, or information
  cannot possibly be copied over, then this function can be used to
  associate that information. One such example is MaxSignal for each of
  the four channels. If this information (MaxSignal) is resident in a
  cmp.h5 file then it will be a column vector containing the maximum
  signal for the channel that was used to call the base. In addition, an
  analyst might wish to see the information from all four channels. This
  function can be used to obtain that information.
}
\seealso{
  \code{\link{PacBioCmpH5-class}}, \code{\link{PacBioPlsH5-class}},
  \code{\link{getBasecallsWithFeatures}}, \code{\link{getAlignmentsWithFeatures}}
}
\examples{
  plsH5 <- PacBioPlsH5(system.file("h5_files", "lambda_example.pls.h5", package = "pbh5"))
  cmpH5 <- PacBioCmpH5(system.file("h5_files", "lambda_example_fromPls.cmp.h5", package = "pbh5"))

  head(associateCmpWithPls(cmpH5, list(plsH5), "QualityValue", idx =
                           c(1:10, 13), collapse = TRUE))

  head(associateCmpWithPls(cmpH5, list(plsH5), "MaxSignal", idx = 1, collapse = TRUE))
}
\keyword{programming}
