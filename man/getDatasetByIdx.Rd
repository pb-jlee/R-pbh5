\name{getDatasetByIdx}
\alias{getDatasetByIdx}
\title{
  Get cmp.h5 Dataset
}
\description{
  Retrieve a dataset by name. Cmp.h5 datasets, such as alignments, are
  stored in a particular format. This function is used to obtain a
  dataset stored in this fashion, see the details section for more
  information.
}
\usage{
  getDatasetByIdx(cmpH5, idx = seq.int(1, nrow(cmpH5)), dsName, convert = NULL)
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to return.
  }
  \item{dsName}{
    The name of the dataset to retrieve.
  }
  \item{convert}{
    A function to apply to each returned slice of the dataset.
  }
}
\details{
  Cmp.h5 datasets are stored in a structured format. For each alignment
  in the alignmentIndex, the columns \code{offsetBegin},
  \code{offsetEnd}, \code{readGroup} define the exact portion of the
  dataset associated with the alignment. This format is described in
  detail in the XXX specification.
}
\value{
  Typically, a list of numeric vectors which are the length of each alignment
  retrieved. In the case of alignments, these are character matrices.
}
\examples{
  cmpH5 <- PacBioCmpH5(system.file("h5_files",
                       "lambda_example_fromPls.cmp.h5", package = "pbh5"))

  z <- getDatasetByIdx(cmpH5, idx = 1:10, dsName = "StartTime", convert =
                        function(a) a - min(a, na.rm = TRUE))
  plot(ecdf(z[[5]]))
}
\keyword{programming}
