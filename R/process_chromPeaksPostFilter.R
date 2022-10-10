#' Extract scan meta and intensity table from XCMSnExp
#' @param x An XCMSnExp object after peak calling
#' @param ms1.noise passed to \link{asIntensityTable}
#' @param ms1.maxPeaks passed to \link{asIntensityTable}
#' @param ms1.maxIdenticalInt passed to \link{asIntensityTable}
#' @param ms2.noise passed to \link{asIntensityTable}
#' @param ms2.maxPeaks passed to \link{asIntensityTable}
#' @param ms2.maxIdenticalInt passed to \link{asIntensityTable}
#' @param ... passed to \code{\link[MSnbase]{spectrapply}} or \code{\link[BiocParallel]{bplapply}}
#' @return a list of "mtab" and "itab"

getMetaIntensityTable <- function(x, 
  ms1.noise = 100, ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 20,
  ms2.noise = 30, ms2.maxPeaks = 100, ms2.maxIdenticalInt = 6, ...
  ) {

  cat("Extracting intensity and scan meta tables ...\n")
  mtable <- asMetaTable(x)
  itable_orig <- asIntensityTable(
    x,  
    ms1.noise = ms1.noise, ms1.maxPeaks = ms1.maxPeaks, ms1.maxIdenticalInt = ms1.maxIdenticalInt,
    ms2.noise = ms2.noise, ms2.maxPeaks = ms2.maxPeaks, ms2.maxIdenticalInt = ms2.maxIdenticalInt, 
    ...) 

  npks <- table(itable_orig$ID)[mtable$ID]
  npks[is.na(npks)] <- 0
  mtable$peakCountFiltered <- npks

  ilist <- split(itable_orig, itable_orig$ID)
  vv <- sapply(seq_len(nrow(mtable)), function(i, itablist) {
    x <- mtable[i, ]
    if (x$msLevel == 2) 
      any(itablist[[x$ID]] < x$precMz - x$precMz * 30 * 1e-6) else
        return(NA)
    }, itablist = ilist)
  mtable$validMS2 <- vv

  list(mtab = mtable, itab = itable_orig)
}

#' post filter of chromPeaks
#' @param x An XCMSnExp object after peak calling
#' @param postfilter a numeric vector of length 2, in for mat c(x, y), meaning at least x peaks with intensity higher than y. 
#' @param mtab scan meta table returned by \link{asMetaTable}
#' @param itab intensity table returned by \link{asIntensityTable}
#' @param BPPARAM arguments passed to bplapply
#' @export
#' 
chromPeaksPostFilter <- function(
  x, mtab, itab, postfilter = c(5, 500), BPPARAM = bpparam()
  ) {

  cat("Chrompeaks post filter ...\n")
  phold <- x@.processHistory 
  peaks <- chromPeaks(x)
  
  fi <- fastmatch::fmatch(itab$ID, mtab$ID)
  iir <- mtab$msLevel[fi] == 1 & itab$intensity >= postfilter[2]
  itable <- itab[iir, ]
  itable$fromFile <- mtab$fromFile[fi][iir]
  itable$rt <- mtab$rt[fi][iir]
  itable$ID <- NULL  
  
  rs <- bplapply(1:nrow(peaks), function(i) {
    i1 <- peaks[i, ]
    mzoff <-i1[["mz"]]*1e-6*25
    irt <- which(dplyr::between(itable$rt, i1[["rtmin"]], i1[["rtmax"]]) & 
                   dplyr::between(itable$mz, i1[["mz"]]-mzoff, i1[["mz"]]+mzoff) & 
                   itable$fromFile == i1[["sample"]])
    ifelse(length(irt) < postfilter[1], FALSE, TRUE)
  }, BPPARAM = BPPARAM)
  rs <- unlist(rs)
  valid <- attr(x, "chromPeakMS2")
  if (!is.null(valid))
    rs <- rs | valid$validMS2
  
  vm <- x@msFeatureData@.xData$chromPeaks[rs, ]  
  rownames(vm) <- xcms:::.featureIDs(nrow(vm), prefix = "CP")
  if (!is.null(valid)) {
    valid <- valid[rs, ]
    rownames(vm) <- rownames(vm)
    attr(x, "chromPeakMS2") <- valid
  }
  chromPeaks(x) <- vm
  x@.processHistory <- phold
  x
}

