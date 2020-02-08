#' extending chrom peaks table with the MS2, chromatogram (with QC)
#' @param x an MSnExp object
#' @param mtab meta table, often returned by \link{asMetaTable}
#' @param itab intensity table, often returned by \link{asIntensityTable}
#' @param ... other parameters passed to bplapply
#' @export
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr between
#' @importFrom xcms chromPeaks


extendedChromPeaks <- function(x, mtab, itab, ...) {
  
  if (!inherits(x, "MSnExp"))
    stop("x should be an object of class MSnExp!")
  
  cpeaks <- xcms::chromPeaks(x)
  cpeaks <- as.data.frame(cpeaks)
  
  cat("Mapping MS2 spectra ...\n")
  mtab2 <- mtab[mtab$msLevel == 2, ]
  
  mtab2_list <- split(mtab2, mtab2$fromFile) 
  l <- bplapply(1:nrow(cpeaks), function(i) {
    c1 <- cpeaks[i, ]
    tmp <- mtab2_list[[as.character(c1[["sample"]])]]
    ic <- which(
      dplyr::between(tmp$rt, c1[["rtmin"]], c1[["rtmax"]] ) & 
        dplyr::between(tmp$precMz, c1[["mzmin"]], c1[["mzmax"]])
    )
    tmp$ID[ic]
  })
  cpeaks$ms2Scan <- sapply(l, paste, collapse = ";")
  
  cat("Evaluating EIC quality ...\n")
  leic <- bplapply(1:nrow(cpeaks), function(i) {
    c1 <- cpeaks[i, ]
    masstol <- c1$mz * 25 * 1e-6
    pk <- xcmsViewerApp::eic(mtab = mtab, itab = itab, file = c1$sample, 
              rt = c(c1$rtmin - 10, c1$rtmax + 10), 
              mz = c(c1$mzmin - masstol, c1$mzmax + masstol))
    if (is.null(pk))
      return(
        structure(rep(NA, 6), names = c("rsq", "rtgap", "intgap", "rtintgap", "truncated", "b"))
      )
    .evalPeak(pk)$stats
  }, ...)
  ltab <- do.call(rbind, leic)
  cpeaks <- cbind(cpeaks, ltab)
  cpeaks$QC <- NA
  cpeaks$QC[cpeaks$rsq > 0.8 & abs(cpeaks$truncated) < 1] <- "+"
  cpeaks
}
