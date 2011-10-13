\name{makeContextDataTable}
\alias{makeContextDataTable}
\alias{associateWithContext}
\alias{summarizeByContext}
\alias{getKmerContext}
\title{
  Compute Statistics Regarding Specific Contexts
}
\description{
  These functions compute statistics like the median IPD for a 
  specific context, or the raw IPD for a single context. 
}
\usage{
  makeContextDataTable(cmpH5, idx = 1:nrow(cmpH5), up = 2, down = 2, fxs
  = list(ipd = getIPD, pw = getPulseWidth, tpos = getTemplatePosition), ...)

associateWithContext(cmpH5, idx = 1:nrow(cmpH5), f = getIPD, up = 2, down = 2, 
    useReference = TRUE, collapse = FALSE)

summarizeByContext(cmpH5, idx = 1:nrow(cmpH5), up = 2, down = 2,
             statF = getIPD, summaryF = function(a) median(a, na.rm =
     T), ...)
}

\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to retrieve.
  }
  \item{up}{
    Number of nucleotides upstream to include in the context.
}
  \item{down}{
    Number of nucleotides downstream to include in the context.
}
  \item{statF}{
    Statistic to be calculated, the default is 'getIPD.'
}
  \item{f}{
    Statistic to be calculated, the default is 'getIPD.'
}
  \item{summaryF}{
    Summary statistic to be calculated, the default is the median.
}
  \item{fxs}{
    A list of statistics to be included in the table, the defaults are
  'getIPD', 'getPulseWidth' and 'getTemplatePosition.'
}
  \item{useReference}{
    Boolean statement whether to use the reference strand.
}  
  \item{collapse}{
    Collapses all context columns into one column.
}
  \item{\dots}{
    Extra arguments.
}
}
\value{
 These functions return data.frame(s) with columns containing either numeric or
 character vectors.
}
\examples{
  ## These functions contain similar information regarding IPD.
  ## This example demonstrates that the median IPD can be computed
     using any one of these functions. 

  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))

  summarize <-summarizeByContext(cmpH5)
  head(summarize)

  associate <- associateWithContext(cmpH5)
  head(associate)
  a5A <- subset(associate, context.P01 == 'A' & context.P02 == 'A' &
                 context.P03 == 'A' & context.P04 == 'A' & context.P05 == 'A')
  head(a5A)
  nrow(a5A)
  median(a5A$elt, na.rm = TRUE)

  makeCDT <- makeContextDataTable(cmpH5)
  m5A <- subset(makeCDT, context.P01 == 'A' & context.P02 == 'A' &
          context.P03 == 'A' & context.P04 == 'A' & context.P05 == 'A')
  median(m5A$elt.ipd, na.rm = TRUE)

  ## Some plots...

  plot(density(log(a5A$elt), na.rm = TRUE), main = "Distribution of IPDs
  for context AAAA")
  boxplot(list(m5A$elt.ipd, m5A$elt.pw), ylim = c(0,1), names = c('IPD','PW'))
}
\keyword{datasets}

