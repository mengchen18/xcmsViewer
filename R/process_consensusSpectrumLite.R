#' calculate consensus spectrum 
#' @param x a data.frame with at least three columns named as "mz", 'intensity' and "ID". 
#'   ID tells which peak is from which spectra.
#' @param ppmtol ppm tolerance
#' @param minProp minimum proportion that a peak should be included
#' @param df whether to return a data.frame, if false, paste into charecter vector
#' 
consensusSpectrumLite <- function(x, ppmtol=20, minProp = 0.3, df = TRUE) {
  nt <- length(unique(x$ID))
  x <- x[order(x$mz), ]
  mz <- x$mz
  mzd <- mz[-1] - mz[-length(mz)]
  i <- which(mzd < pmax(ppmtol*mz[-1]/1e6, 0.003))
  grp_init <- grp <- 1:nrow(x)
  grp[i+1] <- grp_init[i]
  while (any(grp_init != grp)) {
    grp_init <- grp
    grp[i+1] <- grp[i]
  }
  x$grp <- grp
  tg <- table(grp)
  i <- names(tg)[tg >= length(unique(x$ID)) * minProp]
  if (length(i) == 0) {
    if (df)
      r <- data.frame(mz = numeric(0), intensity=numeric(0)) else
        r <- ""
      attr(r, "purity") <- 0
      return(r)
  }
  x <- x[grp %in% i, ]
  pur <- length(unique(x$ID))/nt
  r <- data.frame(
    mz = tapply(x$mz, x$grp, median),
    intensity = tapply(x$intensity, x$grp, sum)
  )
  r$mz <- signif(r$mz, digits = 7)
  if (!df) {
    r <- paste(paste(r$mz, collapse = ";"), paste(r$intensity, collapse = ";"), sep = "|")
  }
  attr(r, "purity") <- pur
  r
}


# #' Estimate PPM tolerance of MS1
# #' @param fd feature data
# #' @return a list of two elements:
# #'  - tol: PPM estimated tolerance
# #'  - fun: weight for a given PPM (based on CDF of PPM)
# #' @importFrom spatstat CDF
# #' 
# .estimateMS1PPMtol <- function(fd, ppmtol = 10) {
  
#   ms1ppmtol <- 1e6*(fd$mzmax - fd$mzmin)/fd$mzmed
#   ms1ppmtol <- ms1ppmtol[ms1ppmtol > 0]
#   ff <- function(x) return(rep(1, length(x)))
  
#   if (length(ms1ppmtol) > 0) {
#     ppmtol <- min(quantile(ms1ppmtol, 0.95), 20)
#     ms1ppmtol <- ms1ppmtol[ms1ppmtol < ppmtol]
#     f <- CDF(density(ms1ppmtol))
#     ff <- function(x) {
#       pmin(rep(1, length(x)), (1/(1-f(1)))*(1-f(x)))
#     }
#   }
#   list(fun = ff, tol = ppmtol)
# }

# #' Function used to weight ms2 score (std lib)
# #' @param x positive integer, the number of fragments
# .sigWeight <- function(x) {
#   if (x <= 0)
#     stop('x needs to be postive integer!')
#   ((-x)/(2*1+abs(x)))^5
# }
