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
  max(0, (observ-expect)/(1-expect))
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
  pp <- .prep_mirrorPlot( peak.upper = measured, peak.lower = standard, ppmtol = ppmtol)
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

#' #' Measure the cosine distance (spectra angle) of two MS2 spectra
#' #' @param scanList  a named list of scans
#' #' @param refList a named list of reference spectra, each list should have three columns:
#' #'   - mz
#' #'   - intensity
#' #'   - IndexMs2
#' #' @param ppmtol ppm tolerance
#' #' 
#' cospecList <- function(scanList, refList, ppmtol) {
#'   
#'   f1 <- scanList$ID
#'   f2 <- refList$IndexMs2
#'   t1 <- scanList$mz
#'   t2 <- refList$mz
#'   ss <- abs(sapply(t2, function(x) (x-t1)/x))
#'   if (!is.matrix(ss))
#'     ss <- matrix(ss, nrow = length(t1))
#'   rr <- which(ss < ppmtol*1e-6, arr.ind = TRUE)
#'   rx <- cbind(f1[rr[, 1]], f2[rr[, 2]])
#'   rx <- unique(rx)
#'   
#'   if (nrow(rx) == 0) {
#'     df <- data.frame(IndexMs2 = character(0), 
#'                      ScanId = character(0),
#'                      Cos = numeric(0))
#'   } else {
#'     x2 <- sapply(1:nrow(rx), function(i) {      
#'       cospec(measured = scanList[scanList$ID == rx[i, 1], ], 
#'              standard = refList[refList$IndexMs2 == rx[i, 2], ], 
#'              ppmtol = ppmtol, alpha = 1)
#'     })
#'     df <- data.frame(
#'       IndexMs2 = rx[, 2],
#'       ScanId = rx[, 1],
#'       Cos = x2,
#'       stringsAsFactors = FALSE
#'     )
#'     df <- df[which(df$Cos > 0), ]
#'   }
#'   attr(df, "scanList") <- unique(scanList$ID)
#'   attr(df, "refList") <- unique(refList$IndexMs2)
#'   df
#' }
#' 
#' 
