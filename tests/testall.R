## Copyright (c) 2010, Pacific Biosciences of California, Inc.

## All rights reserved.

## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:

##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.

##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.

##     * Neither the name of Pacific Biosciences nor the names of its
##       contributors may be used to endorse or promote products derived
##       from this software without specific prior written permission.

## THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED

## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS CONTRIBUTORS
## BE LIABLE FOR ANY

## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
## GOODS OR SERVICES;

## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
## LIABILITY, OR TORT

## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## ###################################################################
##
## Test code for pbh5 classes and API. 
##
## The goals of these tests are to verify the correctness of the
## returned results. One major complicatation is the 0-vs-1 based
## nature of .h5 versus R. The h5r package brings everything into the
## realm of 1-based indexing, but many aspects of the PB .h5 files
## like offsets are encoded as 0-based and therefore much of this
## application level code needs to manage that to present a uniform
## 1-based interface to the user.
##
## ###################################################################
require(pbh5)
source("tinyTestHarness.R")

##
## Make a new TestHarness.
##
TH <- TestHarness()

##
## Compare files
##
cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5", package = "pbh5"))

getDatasetFromCmpH5 <- function(cmpH5, row, ds = "AlnArray") {
  getH5Dataset(cmpH5, paste(cmpH5$alnGroupPath[row], ds, sep = "/"), inMemory = FALSE)
}

getRangeInH5WithOffsets <- function(cmpH5, row, ...) {
  d <- getDatasetFromCmpH5(cmpH5, row, ...)
  d[cmpH5$offsetBegin[row]:cmpH5$offsetEnd[row]]
}

TH("CMP: Alignments Match", all(sapply(1:nrow(cmpH5), function(i) {
  all(getAlignmentsRaw(cmpH5, i)[[1]] == getRangeInH5WithOffsets(cmpH5, i))
})))

TH("CMP: Alignments Match Reverse", {
  alns <- getAlignmentsRaw(cmpH5, nrow(cmpH5):1)
  all(mapply(function(a, b) {
    all(a == b)
  }, alns, sapply(nrow(cmpH5):1, function(i) getRangeInH5WithOffsets(cmpH5, i))))
})

trueAndNA <- function(a, b) {
  all(a[!is.na(a)] == b[!is.na(b)]) && all(is.na(a) == is.na(b))
}

TH("CMP: Alignments And Features 1", {
  alnsAndFeats <- getAlignmentsWithFeatures(cmpH5, features = c("AlnArray", "IPD",
                                                     "PulseWidth", "QualityValue"))
  alns <- getAlignments(cmpH5)
  all(mapply(function(a,b) {
    all(a == b[, c("read", "reference")])
  }, alns, alnsAndFeats))
})

TH("CMP: Alignments And Features 2", {
  alnsAndFeats <- getAlignmentsWithFeatures(cmpH5, features = c("AlnArray", "IPD", "PulseWidth","QualityValue"))
  trueAndNA(getRangeInH5WithOffsets(cmpH5, 1, ds = "IPD"), alnsAndFeats[[1]][,"IPD"])
})

TH("H5: Dataset Error", {
  assertError(getDatasetByIdx(cmpH5, idx = 1:10, dsName = "sdkflsjdfa"))
})

##
## Test the reads in range and coverage.
##
TH("CMP: Test getReadsInRange", {
  require(IRanges)
  
  all(sapply(refGroup(cmpH5)$Name, function(rSeq) {
    N <- 30
    s <- floor(runif(N, 0, 10000))
    e <- floor(s + rexp(N, 1/25))
    aI <- alnIndex(cmpH5)
    aI <- aI[aI$refName == rSeq, ]
    
    all(mapply(function(s, e) {
      x1 <- sort(getReadsInRange(cmpH5, rSeq, s, e))
      x2 <- sort(matchMatrix(findOverlaps(IRanges(s, e), IRanges(aI[, "tStart"], aI[, "tEnd"])))[,"subject"])
      all(cmpH5$ID[x1] == aI$ID[x2])
    }, s, e))
  }))
})

s <- sample(seq.int(1, nrow(cmpH5)), size = nrow(cmpH5)/2)

TH("CMP: Test ipd, pulsewidth",
   (all(sapply(getIPD(cmpH5, idx = s), length) ==  sapply(getPulseWidth(cmpH5, idx = s), length))))

TH("CMP: IPD correct lengths", all(sapply(getIPD(cmpH5), length) == cmpH5$offsetEnd - cmpH5$offsetBegin + 1))

ipdAndPos <- function(cmpH5, idxs) {
  mapply(getTemplatePosition(cmpH5, idxs, F, F),
         getIPD(cmpH5, idxs), FUN = function(a, b) {
           data.frame(pos = a, ipd = b)
         }, SIMPLIFY = FALSE)
}

d <- associateWithContext(cmpH5, 1:10, f = ipdAndPos)
b <- associateWithContext(cmpH5, 1:10)
TH("CMP: IPDs match", all(na.omit(d$elt.ipd == b$elt)))

a <- makeContextDataTable(cmpH5, 1:10)
TH("CMP: context data table match", all(na.omit(a$elt.ipd == b$elt)))

n <- do.call(rbind, getTemplatePosition(cmpH5, withAlignments = TRUE, asDataFrame = TRUE))
n <- subset(n, strand == 0)
TH("CMP: bases match in template position", all(sapply(split(n[,3], n[,1]), function(a) {
  all((sort(table(factor(a, c("A","C","G","T")))) >= 1) == c(F, F, F, T))
})))

ctxts <- getKmerContext(cmpH5, 1:10)
TH("CMP: contexts match", all(sapply(ctxts, length) == sapply(getAlignments(cmpH5, 1:10), nrow)))

blk <- getAlignmentBlock(cmpH5, 1, 20, 29)
TH("CMP: alignment block match", ncol(blk) == 10)

##
## Pulse files
##
basFiles <- list.files(system.file("h5_files", package = "pbh5"), pattern = "bas\\.h5", full.names = TRUE)
basH5s <- lapply(basFiles, PacBioBasH5)
names(basH5s) <- sapply(basH5s, getMovieName)
basH5 <- basH5s[[1]]

TH("PLS basecall lengths", {
  all(sapply(getBasecalls(basH5), length) == getBaseEvents(basH5)[, "numEvent"])
})

x <- doWithPlsAndCmp(cmpH5, basH5s, fx = function(cmpH5, basH5, idxs) {
  mapply(function(bc, alns) {
    list(bc, alns)
  }, getBasecalls(basH5, holeNumbers = cmpH5$holeNumber[idxs]),
         getAlignments(cmpH5, idxs), SIMPLIFY = FALSE)
}, SIMPLIFY = FALSE)



##
## Print and throw any exceptions.
##
TH(action="print")
TH(action="throw")
