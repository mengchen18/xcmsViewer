#' internal function used by mirrorPlotModule
#' @param x data.frame with columne mz and intensity
#' @param col peak colors
#' @param legend.peak.lower the legend for peak.lower intensities
#' @param legend.peak.upper the legend forpeak.upper intensities
#' @param highlight the peak to be peak.lower
#' @param xlim xlim passed to plot
#' @rawNamespace import(graphics, except = layout)
#' @export

.mirrorPlot <- function(x, col, legend.peak.lower="", legend.peak.upper="", highlight=NULL, xlim=NULL) {
  
  if (any(is.infinite(xlim)))
    return()
  if (nrow(x) == 0) 
    return()
  
  if (length(col) == 1)
    col <- rep(col, nrow(x))
  
  if (!is.null(xlim)) {
    ir <- x$mz > xlim[1] & x$mz < xlim[2]
    x <- x[ir, ]
    col <- col[ir]
  }
  
  if (nrow(x) == 0)
    abi <- 1 else
      abi <- max(abs(x$intensity))
    par(mar= c(4, 4, 4, 2))
    plot(x$mz, x$intensity, col = NA, axes = FALSE, xlab = "", ylab = "", ylim = c(-abi, abi), xlim = xlim)
    if (par("usr")[2] - par("usr")[1] < 0.001) {
      if (is.null(xlim))
        mc <- mean(par("usr")[1:2]) else
          mc <- mean(xlim)
        xlim <- c(mc - 0.0005, mc + 0.0005)
        plot(x$mz, x$intensity, col = NA, axes = FALSE, xlab = "", ylab = "", ylim = c(-abi, abi), xlim = xlim)
    }
    
    yax <- par("yaxp")
    abline(h = seq(yax[1], yax[2], length.out = yax[3]+1), lwd = 0.6, col= "gray", lty = 2)
    abline(h = 0, lwd = 0.6)
    ya <- par("yaxp")
    tcks <- seq(ya[1], ya[2], length.out = ya[3]+1)
    tckslab <- tcks
    axis(side = 2, at = tcks, labels = abs(tckslab))
    
    mtext(side = 2, text = "Intensity", line = 3)
    mtext(side = 4, text = "m/z", cex = 0.9, line = 0.5, las = 1)
    
    mtext(side = 1, at = par("usr")[1], legend.peak.lower, adj = 0, cex = 0.9, line = 2)
    mtext(side = 3, at = par("usr")[1], legend.peak.upper, adj = 0, cex = 0.9, line = 2)
    
    if (nrow(x) == 0)
      return()
    
    segments(x$mz, 0, x$mz, x$intensity, col = col)
    if (!is.null(highlight))
      segments(x[highlight, "mz"], 0, x[highlight, "mz"], x[highlight, "intensity"], col = col[highlight], lwd = 3)
    
    i1 <- order(x$intensity)[1:min(3, sum(x$intensity < 0))]
    i2 <- order(x$intensity, decreasing = TRUE)[1:min(3, sum(x$intensity > 0))]
    xlab <- x[c(i1, i2), ]
    
    par(xpd = TRUE)
    text(xlab$mz, xlab$intensity+sign(xlab$intensity)*par("cxy")[2]*0.6, labels = round(xlab$mz, digits = 4), cex = 0.9)
    
}

#' internal function used by mirrorPlotModule
#' @param peak.upper the upper peak
#' @param peak.lower the lower peak
#' @param ppmtol mass tolerance (in ppm) comparing the upper and lower peaks 
.prep_mirrorPlot <- function(peak.upper, peak.lower=NULL, ppmtol = 10) {
  
  if (is.null(peak.lower)) {
    ms <- peak.upper
    cc <- "gray25"
  } else if (is.null(peak.upper)) {
    ms <- peak.lower
    ms$intensity <- -ms$intensity
    cc <- "gray25"
  } else {
    maxm <- max(peak.lower$intensity)
    maxs <- max(peak.upper$intensity)
    mmrat <- maxm/maxs
    
    mdelta <- abs(sapply(peak.lower$mz, function(x) x-peak.upper$mz))
    if (!inherits(mdelta, "matrix"))
      mdelta <- matrix(mdelta, ncol = length(peak.lower$mz))
    ppm <- 1e6 * mdelta / peak.lower$mz    
    ppmtol <- pmax(ppmtol, 1e6*0.002/peak.lower$mz)
    i <- which(ppm < ppmtol, arr.ind = TRUE)
    ppmcut <- ppm[which(ppm < ppmtol)]
    
    peak.upper$mapGroup <-peak.upper$ppm <- NA
    peak.upper$mapGroup[i[, 1]] <- 1:nrow(i)
    peak.upper$ppm[i[, 1]] <- ppmcut
    
    peak.lower$mapGroup <- peak.lower$ppm <- NA
    peak.lower$mapGroup[i[, 2]] <- 1:nrow(i)
    peak.lower$ppm[i[, 2]] <- ppmcut
    
    peak.lower$intensity <- -peak.lower$intensity
    stdn <-peak.upper[, c("mz", "intensity", "mapGroup","ppm")]
    stdn$intensity <- stdn$intensity*mmrat
    ms <- rbind(peak.lower[, c("mz", "intensity", "mapGroup","ppm")], stdn)
    ms <- ms[ms$intensity != 0, ]
    
    cc <- rep("gray25", nrow(ms))
    cc[!is.na(ms$mapGroup)] <- "green"
  }
  list(tab = ms, col = cc)
}


#' Mirror plot of MS2 fragmentation patterns
#' @param peak.lower the peak.lower fragmentation pattern, a data.frame with two columns named "mz" and "intensity"
#' @param peak.upper the peak.upper fragmentation pattern, a data.frame with two columns named "mz" and "intensity"
#' @param legend.peak.lower the legend for peak.lower intensities
#' @param legend.peak.upper the legend forpeak.upper intensities
#' @param ppmtol the mass tolerane used to map the intensity
#' @param highlight the peak to be peak.lower
#' @param xlim xlim passed to plot
#' @export
#' @examples 
#' df <- data.frame(
#'   mz = abs(rnorm(20)*200),
#'   intensity = abs(rnorm(20)*1000)
#' )
#' set.seed(1000)
#' me <- df[sample(1:nrow(df), size = 15), ] + rnorm(5, sd = 0.0005)
#' st <- df[sample(1:nrow(df), size = 10), ]
#' 
#' mirrorPlot(peak.lower = me,
#'            peak.upper = st,
#'            legend.peak.lower  = "xxxxxxxxxxx",
#'            legend.peak.upper = "yyyyyyyyyyyyy")
#' 
#' mirrorPlot(peak.lower = me,
#'            peak.upper = st,
#'            legend.peak.lower  = "xxxxxxxxxxx",
#'            legend.peak.upper = "yyyyyyyyyyyyy",
#'            xlim = c(50, 120))



mirrorPlot <- function(
  peak.upper, peak.lower, 
  legend.peak.lower="", legend.peak.upper="", 
  ppmtol = 10, highlight=NULL, xlim = NULL) {
  if (any(is.infinite(xlim)))
    return(NULL)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  v <- .prep_mirrorPlot(peak.lower = peak.lower, peak.upper = peak.upper, ppmtol)
  .mirrorPlot(v$tab, v$col, legend.peak.lower, legend.peak.upper, highlight = highlight, xlim = xlim )
}

