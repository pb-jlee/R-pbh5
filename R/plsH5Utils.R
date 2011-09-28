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

## XXX: There is a duplication of features which are in both PulseCalls and BaseCalls.
##      you need to deal with that as it makes certain things unaccessible. 

.ofPulses <- c("Channel", "Chi2", "ClassifierQV", "IsPulse", "MaxSignal",
               "MeanSignal", "MidSignal", "StartFrame", "WidthInFrames")

.ofBasecalls <- c("Basecall", "DeletionTag", "DeletionQV",
                  "InsertionQV", "PulseIndex", "PreBaseDeletionQV",
                  "QualityValue", "SubstitutionTag", "SubstitutionQV")

getBasecallsWithFeatures <- function(plsH5, features, group = c( ),
                                     holeNumbers, collapse = FALSE) {
  if (any(features %in% .ofPulses) && class(plsH5) == "PacBioBasH5")
    stop(paste("File:", plsH5@fileName, "doesn't have pulse features."))
  
  if (! all(features %in% c(.ofPulses, .ofBasecalls)))
    stop("Feature must be the correct name of a dataset; possible values:",
         paste(c(.ofPulses, .ofBasecalls), collapse = ", "))
  
  if (missing(holeNumbers))
    holeNumbers <- getHoleNumbers(plsH5)
  
  basecalls <- getBasecall(plsH5, holeNumbers = holeNumbers)
  basecalls <- Filter(Negate(is.null), basecalls)
  holeNumbers <- as.integer(names(basecalls))
  
  if ((hP <- any(features %in% .ofPulses))) {
    pindices  <- getPulseIndex(plsH5, holeNumbers = holeNumbers)
    pfeatures <- intersect(features, .ofPulses)
    
    pdta <- lapply(pfeatures, function(nm) {
      getFromPlsH5(plsH5, "PulseCalls", nm, holeNumbers = holeNumbers)
    })
    
    pres <- lapply(1:length(pindices), function(i) {
      f <- do.call(cbind, lapply(pdta, function(d) {
        if (! is.null(dim(d[[i]]))) {
          d[[i]][pindices[[i]] + 1,, drop=FALSE]
        }
        else {
          d[[i]][pindices[[i]] + 1]
        }
      }))
      colnames(f) <- rep(pfeatures, sapply(pdta, function(d) {
        if (is.null(ncol(d[[1]]))) 1 else ncol(d[[1]])
      }))
      return(f)        
    })
    features <- setdiff(features, pfeatures)
  }

  if (length(features) > 0) 
    res <- lapply(features, function(nm) {
      getFromPlsH5(plsH5, "BaseCalls", nm, holeNumbers = holeNumbers)
    })
  
  res <- lapply(1:length(basecalls), function(i) {
    if (length(features) > 0)
      f <- data.frame(basecalls[[i]], do.call(cbind, lapply(res, function(d) d[[i]])),
                      stringsAsFactors = FALSE)
    else
      f <- data.frame(basecalls[[i]], stringsAsFactors = FALSE)
    
    colnames(f) <- c("basecall", features)
    return(f)
  })

  if (collapse) {
    zmw <- rep(holeNumbers, sapply(res, nrow))
    res <- do.call(rbind, res)
    if (hP)
      res <- cbind(res, do.call(rbind, pres))
    res$zmw <- zmw
  }
  
  else {
    if (hP) {
      res <- mapply(function(a,b) {
        cbind(a,b)
      }, res, pres, SIMPLIFY = FALSE)
    }
    names(res) <- holeNumbers
  }
  return(res)
}


associateCmpWithPls <- function(cmpH5, plsFiles, features,
                                idx = seq.int(1, nrow(cmpH5)), collapse = TRUE) {
  
  associateHole <- function(aln) aln$holeNumber
  
  movieNames <- sapply(plsFiles, function(plsH5) {
    getMovieName(plsH5)
  })

  getPlsFile <- function(nm) {
    plsFiles[[pmatch(nm, movieNames)]]
  }

  fixColnames <- function(a) {
    cn <- table(factor(colnames(a), colnames(a)[!duplicated(colnames(a))]))
    do.call(c, mapply(function(cnt, nm) {
      if (cnt == 1) nm else paste(nm, 1:cnt, sep = ".")
    }, cn, names(cn), SIMPLIFY = FALSE))
  }
  idxByMovie <- split(idx, alnIndex(cmpH5)[idx, "movieName"])
  
  lst <- lapply(idxByMovie, function(idxs) {
    alns <- alnIndex(cmpH5)[idxs,,drop = FALSE]
    plsFile <- getPlsFile(alns[1, "movieName"])
    
    if (length(plsFile) == 0) {
      warning("Wrong pulse files specified.")
      return(NULL)
    }
    featuresForAlns <- getBasecallsWithFeatures(plsFile, features, associateHole(alns))
    alignmentArrays <- getAlignments(cmpH5, idxs)
    
    lapply(1:nrow(alns), function(i) {
      feats  <- featuresForAlns[[i]][alns[i,"rStart"]:alns[i,"rEnd"], -1,
                                     drop = FALSE]
      aArray <- alignmentArrays[[i]]
      c(read = list(aArray[,1]), ref = list(aArray[,2]), 
        lapply(feats, function(a) {
          v <- vector(class(a), nrow(aArray))
          v[aArray[,1] != "-"] <- a
          v[aArray[,1] == "-"] <- NA
          v
        }))
    })
  })

  ## this is necessary to produce a list in the same order as the
  ## input idx.
  lst <- do.call(c, lst)
  names(lst) <- do.call(c, idxByMovie)
  lst <- lst[as.character(idx)]

  if (collapse) {
    oAlnIdx <- rep(as.integer(names(lst)), sapply(lst, function(a) length(a[[1]])))
    someNames <- names(lst[[1]])
    
    lst <- as.data.frame(lapply(1:length(someNames), function(i) {
      do.call(c, lapply(lst, function(l) {
        l[[i]]
      }))
    }), stringsAsFactors = FALSE)
    
    colnames(lst) <- someNames
    lst$alnIdx <- oAlnIdx
    colnames(lst) <- fixColnames(lst)
  }
  
  return(lst)
}

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

  ## based on Pat telling me this was deprecated.
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

setMethod("getSNR", "PacBioPlsH5", function(h5Obj) {
  getH5Dataset(h5Obj, "PulseData/PulseCalls/ZMWMetrics/Snr")[]
})

getNumEvent <- function(plsH5) {
  getH5Dataset(plsH5, "PulseData/PulseCalls/ZMW/NumEvent")[]
}

getReadScore <- function(plsH5) {
  getH5Dataset(plsH5, "PulseData/BaseCalls/ZMWMetrics/ReadScore")[]
}

