#' post filter of chromPeaks
#' @param x An XCMSnExp object after peak calling
#' @param postfilter = c(5, 500), 
#' @param ms1.noise passed to \link{asIntensityTable}
#' @param ms1.maxPeaks passed to \link{asIntensityTable}
#' @param ms1.maxIdenticalInt passed to \link{asIntensityTable}
#' @param ms2.noise passed to \link{asIntensityTable}
#' @param ms2.maxPeaks passed to \link{asIntensityTable}
#' @param ms2.maxIdenticalInt passed to \link{asIntensityTable}
#' @param ... arguments passed bplapply
#' @export
#' 
chromPeaksPostFilter <- function(
  x, postfilter = c(5, 500), 
  ms1.noise = 100, ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 20,
  ms2.noise = 30, ms2.maxPeaks = 100, ms2.maxIdenticalInt = 6, 
  ...) {

  cat("finding chrompeaks ...\n")
  phold <- x@.processHistory 
  peaks <- chromPeaks(x)
  
  cat("Extracting intensity and meta tables ...\n")
  mtable <- asMetaTable(x)
  itable_orig <- asIntensityTable(
    x,  
    ms1.noise = ms1.noise, ms1.maxPeaks = ms1.maxPeaks, ms1.maxIdenticalInt = ms1.maxIdenticalInt,
    ms2.noise = ms2.noise, ms2.maxPeaks = ms2.maxPeaks, ms2.maxIdenticalInt = ms2.maxIdenticalInt, 
    ...)  
  
  fi <- fastmatch::fmatch(itable_orig$ID, mtable$ID)
  iir <- mtable$msLevel[fi] == 1 & itable_orig$intensity >= postfilter[2]
  itable <- itable_orig[iir, ]
  itable$fromFile <- mtable$fromFile[fi][iir]
  itable$rt <- mtable$rt[fi][iir]
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
  x@.processHistory <- phold
  list(mtable = mtable, itable = itable_orig, XCMSnExp = x)
}

