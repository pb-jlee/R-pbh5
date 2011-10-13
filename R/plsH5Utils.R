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

getRegionsTable <- function(plsH5) {
  tbl           <- getH5Dataset(plsH5, "PulseData/Regions")
  regionTypes   <- getH5Attribute(tbl, "RegionTypes")[]
  regionsDescs  <- getH5Attribute(tbl, "RegionDescriptions")[]
  regionSources <- getH5Attribute(tbl, "RegionSources")[]
  columnNames   <- getH5Attribute(tbl, "ColumnNames")[]

  mat <- as.data.frame(tbl[])
  colnames(mat) <- c("holeNumber", "type", "start", "end", "score")

  mat$type  <- regionTypes[mat[,"type"] + 1]
  mat$start <- mat$start + 1
  mat[mat$type != "GlobalAccuracy", ]
}

getAcqParams <- function(plsH5) {
  names <- listH5Contents(g <- getH5Group(plsH5, "ScanData/AcqParams"))$"."$attributes
  s <- sapply(names, function(a) {
    getH5Attribute(g, a)[]
  })
  s$MovieTime <- s$NumFrames / s$FrameRate
  return(s)
}

setGeneric("getNumEvent", function(h5Obj, ...) {
  standardGeneric("getNumEvent")  
})

setMethod("getNumEvent", "PacBioPlsH5", function(h5Obj) {
  getH5Dataset(h5Obj, "PulseData/PulseCalls/ZMW/NumEvent")[]
})

setMethod("getNumEvent", "PacBioBasH5", function(h5Obj) {
  getH5Dataset(h5Obj, "PulseData/BasCalls/ZMW/NumEvent")[]
})

getReadScore <- function(basH5) {
  getH5Dataset(basH5, "PulseData/BaseCalls/ZMWMetrics/ReadScore")[]
}

setMethod("getSNR", "PacBioBasH5", function(h5Obj) {
  snr <- getH5Dataset(h5Obj, "PulseData/BaseCalls/ZMWMetrics/HQRegionSNR")[]
  colnames(snr) <- getChannelToBaseMap(h5Obj)
  snr
})

setMethod("getSNR", "PacBioPlsH5", function(h5Obj) {
  snr <- getH5Dataset(h5Obj, "PulseData/PulseCalls/ZMWMetrics/Snr")[]
  colnames(snr) <- getChannelToBaseMap(h5Obj)
  snr
})
