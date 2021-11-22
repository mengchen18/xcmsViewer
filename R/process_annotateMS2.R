#' Annotate using MS2 spectra
#' @param x a prunedXcmsSet object
#' @param ms1Annot the MS1 annotation table returned by \code{annotateMS1}
#' @param ppmtol the mass tolerence given by parts per million (PPM)
#' @import parallel
#' @import BiocParallel
#' @import fastmatch

annotateMS2 <- function(x, ms1Annot, ppmtol = 20) {
  
  fdata <- fData(x@featureSet)
  im <- which(ms1Annot$ms2_purity > 0.5)
  
  v <- sapply(im, function(i) {
    cos <- NA
    fid <- ms1Annot$ID[i]
    measured_cs <- fdata$ms2spectrum[fdata$ID == fid]
    if (nchar(measured_cs) > 0) {
      measured_cs <- str2spectra(measured_cs)
      # measured_cs <- strsplit(measured_cs, split = "\\|")[[1]]
      # measured_cs <- lapply(measured_cs, function(it) strsplit(it, ";")[[1]])
      # measured_cs <- data.frame(
      #   mz = as.numeric(measured_cs[[1]]),
      #   intensity = as.numeric(measured_cs[[2]]))
      if (nchar(ms1Annot$ms2_mass[i]) > 0) {
        ref_cs <- data.frame(
          mz = as.numeric(strsplit(ms1Annot$ms2_mass[i], ";")[[1]]),
          intensity = as.numeric(strsplit(ms1Annot$ms2_intensity[i], ";")[[1]])
        )
        cos <- cospec(measured = measured_cs, standard = ref_cs)
      }
    }
    attr(cos, "ppmtol") <- ppmtol
    cos
  })
  
  ms1Annot$ms2_cos <- NA
  ms1Annot$ms2_cos[im] <- v
  ms1Annot
}


scoreAnnot <- function(x, an2) {
  fdata <- fData(x@featureSet)
  ### scoring annotation
  # MS2 score
  nmea <- stringr::str_count(fdata$ms2spectrum, pattern = ";")/2 + as.integer(nchar(fdata$ms2spectrum) > 0)
  names(nmea) <- rownames(fdata)
  nref <- stringr::str_count(an2$ms2_mass, pattern = ";") + as.integer(nchar(an2$ms2_mass) > 0)
  nref[is.na(an2$ms2_intensity)] <- NA
  scr_ms2 <- ms2Score(n_measured=nmea[an2$ID], n_ref = nref, cos = an2$ms2_cos, std = an2$internalStd)
  
  # MS1 score
  scr_ms1 <- ms1Score(x, deltaPPM = an2$DeltaPPM)
  
  # RT score
  rt <- structure(fdata$rtmed, names = rownames(fdata))
  rtstd <- an2$RT
  rtstd[!an2$internalStd] <- NA
  scr_rt <- rtScore(rt[an2$ID], an2$RT)
  
  # Final score
  an2$score_ms1 <- scr_ms1
  an2$score_ms2 <- scr_ms2
  an2$score_rt <- scr_rt
  an2$Score <- (scr_ms1 + scr_ms2 + scr_rt) * an2$IPS + an2$primary/2
  an2
}

# MS2 score of internal/public
# MS1 score 
# RT score (only internal)
# primary bonus * 1 / extended * 0.75

rtScore <- function(rt, rtStd, maxScore = 60*15) {
  scr <- rep(0, length(rtStd))
  nd <- abs(rt - rtStd)/maxScore
  scr[which(nd > 0.25)] <- -0.25
  scr[which(nd < 0.1)] <- 0.75
  scr[which(nd < 0.05)] <- 1
  scr
}

ms2Score <- function(n_measured, n_ref, cos, std) {
  scr_ms2 <- cos
  scr_ms2[is.na(scr_ms2)] <- 0
  
  i <- which(std & n_measured == 0 & n_ref == 0)
  scr_ms2[i] <- 0.25
  i <- which(std & ((n_measured >= 2 & n_ref == 0) | (n_measured == 0 & n_ref >= 2)))
  scr_ms2[i] <- scr_ms2[i] - 0.25
  i <- which(std & ((n_measured >= 6 & n_ref == 0) | (n_measured == 0 & n_ref >= 6)))
  scr_ms2[i] <- scr_ms2[i] - 0.5
  scr_ms2
}

#' @importFrom stats mad pnorm quantile
ms1Score <- function(x, deltaPPM, lower = 5, upper = 10) {
  fdata <- fData(x@featureSet)
  sds <- sapply(fdata$peakidx, function(i) {
    mz <- x@peak@table$mz[i]
    mad(1e6*mz/median(mz))
  })
  sds[sds==0] <- NA
  sds <- quantile(sds, na.rm = TRUE, probs = 0.75)*2
  sds <- min(upper, sds)
  sds <- max(lower, sds)
  pp <- pnorm(q = deltaPPM, sd = sds, mean = 0, lower.tail = FALSE)
  2*pp
}
