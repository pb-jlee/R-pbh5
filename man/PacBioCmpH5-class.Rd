\name{PacBioCmpH5-class}
\Rdversion{1.1}
\docType{class}
\alias{PacBioCmpH5-class}
\alias{H5ObjOrNull-class}
\alias{matrixOrNull-class}
\alias{stringOrNull-class}
\alias{head,PacBioCmpH5-method}
\alias{nrow,PacBioCmpH5-method}
\alias{show,PacBioCmpH5-method}
\alias{summary,PacBioCmpH5-method}
\alias{getMovieName,PacBioCmpH5-method}
\alias{$,PacBioCmpH5-method}
\alias{[,PacBioCmpH5,ANY,ANY,ANY-method}
\alias{PacBioCmpH5}
\alias{alnGroup}
\alias{refGroup}
\alias{refInfo}
\alias{movieInfo}
\alias{getSNR,PacBioCmpH5-method}
\alias{getEviconsCalls}
\alias{getEviconsConfidence}

\title{Class "PacBioCmpH5"}
\description{
  The PacBioCmpH5 (pronounced PacBio Comp H5) class represents
  alignments to reference sequences. Additionally, pulse information as
  well as quality values can be stored within this file.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{PacBioCmpH5(fileName)}.
}
\section{Slots}{
  \describe{
    \item{\code{AlnIndex}:}{Object of class \code{"data.frame"}
      The alignment index represents each alignment with a single
      row. The row gives information about the position, quality, and
      strand of the alignment. Additionally, information about how
      to access read/ref bases is stored in the row. 
    }
    \item{\code{AlnGroup}:}{Object of class \code{"data.frame"}
      representing the alignment groups in the file. An alignment group
      can be used to partition the reads into different categories in a
      structured fashion, e.g., machines or movies. Typically, this
      represents reads coming from the same movie.
    }
    \item{\code{RefGroup}:}{Object of class \code{"data.frame"}
      representing the reference sequences with alignments in the
      file. This object contains information about which reads map to
      which references as well as some information about the reference
      sequence itself.
    }
    \item{\code{MovieInfo}:}{Object of class \code{"data.frame"}
      representing information about the movies used during the alignment.
    }
    \item{\code{RefInfo}:}{Object of class \code{"data.frame"}
      representing information about the references used during the
      alignment. This data.frame will contain all references which were
      used in the alignment process, whereas the RefGroup contains only
      those references which had one or more reads with an alignment.
    }
    \item{\code{isSorted}:}{Object of class logical}
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
    \item{head}{\code{signature(x = "PacBioCmpH5")}: ... }
    \item{nrow}{\code{signature(x = "PacBioCmpH5")}: ... }
    \item{show}{\code{signature(object = "PacBioCmpH5")}: ... }
  }
}
\examples{
## open a handle to a cmp.h5 file.
cmpH5 <- PacBioCmpH5(system.file("h5_files", "aligned_reads.cmp.h5", package = "pbh5"))

## print a short description of the data.
show(cmpH5)

## get the contents of the entire file.
contents <- listH5Contents(cmpH5)

## contents is a list with an element for each object (dataset, group)
names(contents)

## the alignment index is the core of the file. The alignment index
## contains information about each alignment. This information can be
## used to compute summary statistics on the data.
head(alnIndex(cmpH5))

## direct access to the alignmentIndex is usually not necessary as
## most a large number of functions are available for accessing data
## associated with an alignment, e.g.,
plot(density(getAccuracy(cmpH5)))

## plot density of read length
plot(density(getReadLength(cmpH5)), log = 'x')

##
## coverage plots
##
## coverage plots take a reference sequence name. These are sanitized
## strings representing the reference sequence. These will typically
## be chromosomes, but may be otherwise.
##
cvg <- getCoverageInRange(cmpH5,  1)

## the coverage vector will be the SequenceLength.
stopifnot(length(cvg) ==  refGroup(cmpH5)$Length)

## summarize coverage vector.
summary(cvg)

## plot coverage vector
plot(cvg, type = 'l', col = 'gray')
lines(supsmu(1:length(cvg), cvg, span = .01), col = 'red')

## retrieve all of the alignments
alns <- getAlignments(cmpH5)
head(alns[[1]])

## Compute mismatch-insertion-deletion
a <- do.call(rbind, alns)
mosaicplot(prop.table(table(read = a[,1], reference = a[,2])))

## Some of the most useful functionality is via the
## getAlignmentsWithFeatures function.
aAndF <- getAlignmentsWithFeatures(cmpH5, features = c("IPD"), idx = 1, collapse = TRUE)
head(aAndF)
boxplot(IPD ~ reference, data = aAndF, log = 'y', ylim = c(.05, 10))
}
\keyword{classes}
