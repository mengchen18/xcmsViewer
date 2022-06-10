#' Filter MS spectrum
#' @description Filtering noise in a spectrum. 
#' @param x An Spectrum object
#' @param noise the noise level
#' @param maxPeaks maximum number of peaks to retain
#' @param maxIdenticalInt maximum number of idential intensities (often noise)
#' @param maxMZ maximum MZ, peaks at higher mz will be removed first, used when filter MS2 
#'  spectra to limit mz lower than mass of parent ion
#' @export

filterSpectrum <- function(x, noise = 30, maxPeaks = 100, maxIdenticalInt = 5, maxMZ = Inf) {
  
  if (!inherits(x, "Spectrum"))
    stop("x should be an object of Spectrum object!")
  
  imz <- which(x@mz < maxMZ)
  x@mz <- x@mz[imz]
  x@intensity <- x@intensity[imz]
  
  it <- table(x@intensity)
  it <- which(it >= maxIdenticalInt)
  if (length(it) != 0)
    noise <- max(noise, max(as.integer(names(it)))) else
      noise <- noise
  i <- which(x@intensity > noise)
  np <- length(i)
  
  if (np > maxPeaks) {
    oi <- order(x@intensity[i], decreasing = TRUE)
    oi <- oi[-(1:maxPeaks)]
    i <- i[-oi]
  }
  
  data.frame(
    mz = x@mz[i],
    intensity = x@intensity[i]
  )
}



#' summarize MSnExp into intensity list, intensities are filtered
#' @param x an MSnExp object
#' @param ms1.noise noise level of ms1 spectra, passed to \code{\link{filterSpectrum}}
#' @param ms1.maxPeaks maximum number of retained peaks of ms1 spectra, passed to \code{\link{filterSpectrum}}
#' @param ms1.maxIdenticalInt maximum number of identical intensities allowed of ms1 spectra, passed to \code{\link{filterSpectrum}}
#' @param ms2.noise noise level of ms2 spectra, passed to \code{\link{filterSpectrum}}
#' @param ms2.maxPeaks maximum number of retained peaks of ms2 spectra, passed to \code{\link{filterSpectrum}}
#' @param ms2.maxIdenticalInt maximum number of identical intensities allowed of ms2 spectra, passed to \code{\link{filterSpectrum}}
#' @param ... passed to \code{\link[MSnbase]{spectrapply}} or \code{\link[BiocParallel]{bplapply}}
#' @importFrom MSnbase spectrapply
#' @importFrom BiocParallel bplapply
#' @export

asIntensityTable <- function(
  x, ms1.noise = 100, ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 20,
  ms2.noise = 30, ms2.maxPeaks = 100, ms2.maxIdenticalInt = 6, ...
) {
  
  if (!inherits(x, "MSnExp"))
    stop("x should be an object of class MSnExp!")
  
  on.exit(gc())
  
  xl <- spectrapply(x, function(xx) {
    mslv <- MSnbase::msLevel(xx)
    if (mslv == 1) {
      noise <- ms1.noise
      maxPeaks <- ms1.maxPeaks
      maxIdenticalInt <- ms1.maxIdenticalInt
      maxMz <- Inf
    } else {
      noise <- ms2.noise
      maxPeaks <- ms2.maxPeaks
      maxIdenticalInt <- ms2.maxIdenticalInt    
      maxMz <- xx@precursorMz*3
    }
    filterSpectrum(
      xx, noise = noise, maxPeaks = maxPeaks, maxIdenticalInt = maxIdenticalInt, maxMZ=maxMz
    )
  }, ...)
  df <- do.call(rbind, xl)
  df$ID <- gsub("\\.[0-9]+$", "", rownames(df))
  df
}

#' quick summary of MSnExp
#' @param x an MSnExp object or an oject inheriting it.
#' @export
#' 
asMetaTable <- function(x) {
  
  if (!inherits(x, "MSnExp"))
    stop("x should be an object of class MSnExp!")
  
  on.exit(gc())
  
  .d <- x@featureData@data
  df <- data.frame(
    scanNum = .d$spIdx,
    acquisitionNum = .d$acquisitionNum,
    rt = .d$retentionTime,
    tic = .d$totIonCurrent,
    peakCount = .d$originalPeaksCount,
    msLevel = .d$msLevel,
    fromFile = .d$fileIdx,
    precScanNum = .d$precursorScanNum,
    precMz = .d$precursorMZ,
    precCharge = .d$precursorCharge,
    precIntensity = .d$precursorIntensity,
    ID = rownames(.d),
    stringsAsFactors = FALSE
  )
  i <- which(df$msLevel == 1)
  df$precScanNum[i] <- NA
  df$precMz[i] <- NA
  df$precCharge[i] <- NA
  df$precIntensity[i] <- NA
  df
}