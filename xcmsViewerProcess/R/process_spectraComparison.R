#' normalized cosine distance between two positve vectors
#' @param a a non-negative numerical vector
#' @param b the second non-negative numerical vector has the same length as a
#' @param n power factor for a and b, a and b will be transformed to a^n and b^n
#' @details 
#'  if a and b are non-negative, the expected cosine distance between a and b will
#'  not be 0. Assuming the elements in a and b are random, the normalized cosine
#'  distance is calculated as:
#'    observ <- sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2)))
#'    expect <- sum(mean(a)*mean(b)*length(a))/(sqrt(sum(a^2))*sqrt(sum(b^2)))
#'    (observ-expect)/(1-expect)
#' @return the normalized cosine distance
#' @examples 
#'  x <- abs(rnorm(50))
#'  y <- abs(rnorm(50))
#'  cosnorm(x, y)
#' @export

cosnorm <- function(a, b, n=1) {
  la <- length(a)
  if (any(a < 0) || any(b < 0) || 
      la != length(b) || 
      !is.numeric(a) || !is.numeric(b))
    stop("a and b should non-negative numerical vectors with the same length")
  
  a <- a^n
  b <- b^n
  dn <- sqrt(sum(a^2))*sqrt(sum(b^2))
  observ <- sum(a*b)/dn
  expect <- mean(a)*mean(b)*la/dn
  (observ-expect)/(1-expect)
}


#' Measure the cosine distance (spectra angle) of two MS2 spectra
#' @param measured the measured spectra
#' @param standard the standard spectra compare against
#' @param ppmtol ppm tolerance
#' @param alpha a numerical value. the power applied to the peak intensity, usually between 0-1 to give less weight 
#'  on super high intensity peak (therefore, more weights on the number of matched peaks)
#' @importFrom reshape2 acast
#' @export

cospec <- function(measured, standard, ppmtol = 10, alpha = 1) {
  pp <- xcmsViewerApp::.prep_mirrorPlot( peak.upper = measured, peak.lower = standard, ppmtol = ppmtol)
  tab <- pp$tab
  tab$sig <- sign(tab$intensity)
  i <- !is.na(tab$mapGroup)
  tab$mz[i] <- tab$mapGroup[i]
  mz <- as.character(tab$mz)
  umz <- unique(mz)
  ss <- structure(tab$intensity, names = mz)
  mm <- rbind(
     ss[tab$sig == -1][umz],
     ss[tab$sig == 1][umz]  
   )
  mm[is.na(mm)] <- 0
  mm <- abs(mm)^alpha
  cosnorm(mm[1, ], mm[2, ])
}


#' mapping MS2 spectra
#' @param ftid feature id to be explored
#' @param possibleRefs possible references, return by annotation by mass (querymass)
#' @param refSpectra reference spectra, a list has two elements named as 'meta' and 'peakList',e.g. something returned by "parseHMDBMSMSXML"
#' @param peakidx the peak idx of that ftid
#' @param peakTab the peak table
#' @param itab the itab
#' @param mode either "neg" or "pos" for negative and positive mode
#' @importFrom fastmatch %fin%
#' @importFrom fastmatch fmatch

spectraCosDist <- function(
  ftid,
  possibleRefs, 
  refSpectra, 
  peakidx,
  peakTab,
  itab,
  mode = c("pos", "neg")[1]) {
  
  mode <- switch(
    mode, 
    pos = "positive",
    neg = "negative")
  
  possibleRefSpectraIDs <- possibleRefs[[ftid]]$accession
  if (is.null(possibleRefSpectraIDs))
    return(NULL)
  
  meta <- refSpectra$meta
  meta$`ionization-mode` <- tolower(meta$`ionization-mode`)
  i <- meta$`ionization-mode` %fin% c(mode, "n/a") & meta$`database-id` %fin% possibleRefSpectraIDs
  
  if (!any(i))
    return(NULL)
  
  meta_s <- meta[i, ]
  an_peaks <- refSpectra$peakList[meta_s$id2]
  
  # peakidx <- featureTab[ftid, "peakidx"][[1]]
  scanids <- unlist(strsplit(peakTab[peakidx, "ms2Scan"], ";"))
  if (length(scanids) == 0)
    return(NULL)
  qe_peaks <- itab[itab$ID %fin% scanids, ]
  qe_peaks <- split(qe_peaks, qe_peaks$ID)
  
  cm <- sapply(qe_peaks, function(x) {
    sapply(an_peaks, function(y) cospec(x, y, alpha = 0.5))
  })
  if (!is.matrix(cm))
    cm <- matrix(
      cm, nrow = length(an_peaks), 
      dimnames = list(names(an_peaks), names(qe_peaks))
      )

  
  ii <- which(cm > 0.3, arr.ind = TRUE)
  if (length(ii) == 0)
    return(NULL)

  df <- data.frame(
    annot_peaks = rownames(cm)[ii[, 1]],
    query_peaks = colnames(cm)[ii[, 2]],
    cos = cm[ii],
    stringsAsFactors = FALSE)
  df <- df[order(df$cos, decreasing = TRUE), ]
  df$"database-id" <- meta_s[fmatch(df$annot_peaks, meta_s$id2), "database-id"]
  df
}
