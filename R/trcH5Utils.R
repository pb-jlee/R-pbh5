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

doWithPlsAndTrc <- function(plsH5, trcH5, holeNumbers, fx = rowMeans) {
  holes <- getH5Dataset(trcH5, "TraceData/HoleNumber")[]
  mHoles <- match(holeNumbers, holes)
  trcs <- getTraces(trcH5)[mHoles, ]

  mapply(function(sf, wInF, hole) {
    trc <- trcs[hole,,]

    do.call(rbind, mapply(function(s, w) {
      fx(trc[,s:w])
    }, sf, wInF, SIMPLIFY = FALSE))
    
  }, getStartFrame(plsH5, holeNumbers = holeNumbers),
         getWidthInFrames(plsH5, holeNumbers = holeNumbers), as.list(holeNumbers),
         SIMPLIFY = FALSE)
}

makePlsAndTrcProcessor <- function(plsH5, trcH5, fx = rowMeans, ..., .holeNumbers = NULL, .nBlocks = 50) {
  holes   <- if (is.null(.holeNumbers)) getH5Dataset(trcH5, "TraceData/HoleNumber")[] else .holeNumbers
  data    <- NULL
  buckets <- split(holes, cut(holes, .nBlocks))
  
  i <- 1
  j <- 1
  
  function() {
    if (i > length(buckets)) return(NULL)
    if (is.null(data)) {
      data <<- doWithPlsAndTrc(plsH5, trcH5, holeNumbers = buckets[[i]], fx = fx)
      i    <<- i + 1
      j    <<- 1
    }
    res <- data[[j]]
    j   <<- j + 1

    if (j > length(data)) {
      data <<- NULL
    }

    return(res)
  }
}

applyZMWs <- function(trcH5, fx, ..., .whichZMWs = NULL, .blockSize = 100) {
  trcs <- getTraces(trcH5)

  if (! missing(.whichZMWs)) {
    nZMWs <- length(.whichZMWs)
    wZMWs <- sort(.whichZMWs)
    oZMWs <- order(.whichZMWs)
  } else {
    nZMWs <- dim(trcs)[1]
    wZMWs <- seq.int(1, nZMWs)
  }
  
  nBlocks <- round(nZMWs/.blockSize)
  sV <- rep(1:nBlocks, each = .blockSize)
  
  if (length(sV) > nZMWs) {
    sV <- sV[1:nZMWs]
  } else if (length(sV) < nZMWs) {
    sV <- c(sV, rep(sV[length(sV)], nZMWs - length(sV)))
  }
  stopifnot(length(sV) == nZMWs)
  
  stats <- do.call(c, tapply(wZMWs, sV, function(blck) {
    zmws <- trcs[blck,,drop=FALSE]
    lapply(1:(dim(zmws)[1]), function(i) fx(zmws[i,,], ...))
  }))
  
  if (! missing(.whichZMWs)) {
    names(stats) <- .whichZMWs[oZMWs]
  } else {
    names(stats) <- wZMWs
  }
  stats <- stats[as.character(wZMWs)]
  return(stats)
}
