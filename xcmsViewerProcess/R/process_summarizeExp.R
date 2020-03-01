#' Summarize XCMSnExp data
#' @param x an object of XCMSnExp
#' @param mode mode. pos or neg
#' @param massTab a data.frame of mass table having at least one column named as "monoisotopic_molecular_weight"
#' @param refSpectra reference MS2 spectra used annotate the experimental MS2 spectra
#' @param ms1.noise passed to \link{asIntensityTable}
#' @param ms1.maxPeaks passed to \link{asIntensityTable}
#' @param ms1.maxIdenticalInt passed to \link{asIntensityTable}
#' @param ms2.noise passed to \link{asIntensityTable}
#' @param ms2.maxPeaks passed to \link{asIntensityTable}
#' @param ms2.maxIdenticalInt passed to \link{asIntensityTable}
#' @param QC whether the QC of EIC should be calculated, passed to \link{extendedChromPeaks}
#' @param ... other parameters passed to bplapply 
#' @importFrom fastmatch %fin%
#' @import MAIT
#' @import xcms
#' @import BiocParallel
#' @export
#' 
summarizeExp <- function(x, mode = c("pos", "neg")[1], massTab=NULL, refSpectra = NULL, 
  ms1.noise = 100, ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 20,
  ms2.noise = 30, ms2.maxPeaks = 100, ms2.maxIdenticalInt = 6, QC = TRUE, ...) {

  xl <- xcms_summarize(x = x, 
    ms1.noise = ms1.noise, ms1.maxPeaks = ms1.maxPeaks, ms1.maxIdenticalInt = ms1.maxIdenticalInt,
    ms2.noise = ms2.noise, ms2.maxPeaks = ms2.maxPeaks, ms2.maxIdenticalInt = ms2.maxIdenticalInt, 
    QC = QC, ...)
  
  xcms_annotate(xl, mode = mode, massTab=massTab, refSpectra = refSpectra, ...)
}