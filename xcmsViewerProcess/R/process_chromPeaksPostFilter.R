#' post filter of chromPeaks
#' @param x An XCMSnExp object afater peak calling
#' @param postfilter = c(5, 500), 
# ' @param ms1.noise
# ' @param ms1.maxPeaks
# ' @param x, ms1.noise = 100, 
# ' @param ms1.maxPeaks = Inf, 
# ' @param ms1.maxIdenticalInt = 20,
# ' @param ms2.noise = 30, 
# ' @param ms2.maxPeaks = 100, 
# ' @param ms2.maxIdenticalInt = 6,
#' @param ... arguments passed bplapply
#' @export
#' 
chromPeaksPostFilter <- function(x, postfilter = c(5, 500), ...) {
  cat("finding chrompeaks ...\n")
  
  peaks <- chromPeaks(x)
  
  cat("Extracting intensity and meta tables ...\n")
  mtable <- asMetaTable(x)
  itable <- asIntensityTable(
    x, ms1.noise = postfilter[2], ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 100, ms2.noise = Inf, ...
  )
  fi <- fastmatch::fmatch(itable$ID, mtable$ID)
  itable$fromFile <- mtable$fromFile[fi]
  itable$rt <- mtable$rt[fi]
  itable$ID <- NULL
  
  cat("Post filtering ...\n")
  rs <- bplapply(1:nrow(peaks), function(i) {
    i1 <- peaks[i, ]
    mzoff <-i1[["mz"]]*1e-6*25
    irt <- which(dplyr::between(itable$rt, i1[["rtmin"]], i1[["rtmax"]]) & 
                   dplyr::between(itable$mz, i1[["mz"]]-mzoff, i1[["mz"]]+mzoff) & 
                   itable$fromFile == i1[["sample"]])
    ifelse(length(irt) < postfilter[1], FALSE, TRUE)
  }, ...)
  rs <- unlist(rs)
  
  vm <- x@msFeatureData@.xData$chromPeaks[rs, ]
  rownames(vm) <- xcms:::.featureIDs(nrow(vm), prefix = "CP")
  chromPeaks(x) <- vm
  x
}
