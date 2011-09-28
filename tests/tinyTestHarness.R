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


## #########################################################################
##
## A simple test harness. See testall.R for usage.
##
##
## #########################################################################

assertError <- function(expr) {
  tryCatch({{expr}; FALSE}, simpleError = function(e) {
    return(TRUE)
  })
}

TestHarness <- function() {
  tests <- list()

  getTime <- function(elt) {
    elt[["time"]][3]
  }
  getResult <- function(elt) {
    elt[["result"]]
  }
  printResults <- function() {
    mwidth <- max(nchar(names(tests))) + 5
    fmtString <- paste("\t%-", mwidth, "s %-10g %-10s\n", sep = "")

    cat(sprintf("%s Results for %d tests %s \n\n", paste(rep("-", 30), collapse = ""), length(tests),
                paste(rep("-", 30), collapse = "")))
    
    for (elt in names(tests)) {
      cat(sprintf(fmtString, elt, getTime(tests[[elt]]), getResult(tests[[elt]])))
    }
  }
  
  function(nm, test, action = c("test", "print", "throw")) {
    action <- match.arg(action)
    switch(action,
           test = {
             tm <- system.time({
               b <- tryCatch(test, simpleError = function(e) {
                 return(FALSE)
               })
             })
             tests[[nm]] <<- list("result" = b, "time" = tm)
           },
           print = {
             printResults()
           },
           throw = {
             errs <- ! sapply(tests, getResult)
             if (any(errs)) {
               stop(simpleError(paste("Tests in error:\n",
                                      paste(paste("\t", names(tests)[errs], sep = ""),
                                            collapse = "\n"), sep = "")))
             }
           })
  }
}
