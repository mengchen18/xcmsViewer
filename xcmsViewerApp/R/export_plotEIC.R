#' prepare input for plotEIC function
#' @param x a simple intensity table representing chromatogram, it should be 
#'   a data.frame with three columns named "mz", "intensity" and file (integer)
#' @param diffPheno a vector to indicate which file belongs to which phenotype
#' @param select integer; which file is selected, will be shown in red
#' @export
#' @examples 
#' x <- readRDS(system.file(package="xcmsViewerApp", "extdata", "00_exampleData_processed.RDS"))
#' mtable <- x$scanMetaTab
#' itable <- x$scanIntensityTab
#' chrom <- eic(itab=itable, mtab=mtable, rt = c(0, Inf), mz=c(225, 255), na.rm = TRUE)
#'
#' cc1 <- colCode(chrom, diffPheno=c(1:2), select = NULL) 
#' plotEIC(chrom, unifile=cc1$unifile, col=cc1$col, diffPheno=cc1$diffPheno, vline=NULL)
#'
#' cc2 <- colCode(chrom, diffPheno=c(1:2), select = 1)   
#' plotEIC(chrom, unifile=cc2$unifile, col=cc2$col, diffPheno=cc2$diffPheno, vline=NULL)

colCode <- function(x, diffPheno=NULL, select = NULL) {

  unifile <- unique(x$file)
  uniquePheno <- unique(diffPheno)
  
  if (is.null(select)) {
    if (!is.null(diffPheno) && length(diffPheno) > 0) {
      set.seed(13)
      col <- randomcoloR::distinctColorPalette(length(uniquePheno))
      names(col) <- uniquePheno
    } else {
      col <- c("EIC"="black")
      diffPheno <- rep("EIC", length(unifile))
    }
  } else {
    # if (!is.null(diffPheno))
    #   message("Argument 'selected' is given, diffPheno is ignored.")
    col <- c(unselected = "gray75", selected = "red")
    diffPheno <- rep("unselected", length(unifile))
    diffPheno[unifile %fin% select] <- "selected"
  }

  list(unifile = unifile, col = col, diffPheno = diffPheno)
}


#' Plot extracted ion chromatogram
#' @param x a simple intensity table representing chromatogram, it should be 
#'   a data.frame with three columns named "mz", "intensity" and file (integer)
#' @param unifile the unique files
#' @param col color for each line
#' @param diffPheno a vector to indicate which file belongs to which phenotype
#' @param vline a numerical vector to specify vertical lines to draw. If the name of 
#' of the vector is "rtmin", "rtmax", the two line will be omitted, instead a green 
#' area will be drawn.
#' @param ... other parameter passed to plot
#' @export
#' @examples 
#' x <- readRDS(system.file(package="xcmsViewerApp", "extdata", "00_exampleData_processed.RDS"))
#' mtable <- x$scanMetaTab
#' itable <- x$scanIntensityTab
#' chrom <- eic(itab=itable, mtab=mtable, rt = c(0, Inf), mz=c(225, 255), na.rm = TRUE)
#'
#' cc1 <- colCode(chrom, diffPheno=c(1:2), select = NULL) 
#' plotEIC(chrom, unifile=cc1$unifile, col=cc1$col, diffPheno=cc1$diffPheno, vline=NULL)
#'
#' cc2 <- colCode(chrom, diffPheno=c(1:2), select = 1)   
#' plotEIC(chrom, unifile=cc2$unifile, col=cc2$col, diffPheno=cc2$diffPheno, vline=NULL)
#' 
plotEIC <- function(x, unifile, col, diffPheno, vline=NULL, ...) {
  
  par(mar = c(4, 4, 1, 1))

  plot(x$rt, x$intensity, col = NA, ylab = "Intensity", xlab = "Retention time (second)", axes = FALSE, ...)
  
  if (!is.null(vline)) {
    rect(min(vline), min(x$intensity), max(vline), max(x$intensity), 
      col=rgb(152, 251, 152, maxColorValue = 255, alpha = 25), border = NA)
    abline(v = vline[setdiff(names(vline), c("rtmin", "rtmax"))], lty = 2)
  }
  axis(1)
  axis(2)
  yaxp <- par("yaxp")
  abline(h = seq(yaxp[1], yaxp[2], length.out = yaxp[3]), lty = 3, col = "gray75")
  for (i in unifile[diffPheno != "selected"]) {
    ii <- x$file == i
    lines(x$rt[ii], x$intensity[ii], col = col[diffPheno[match(i, unifile)]])
  }
  for (i in unifile[diffPheno == "selected"]) {
    ii <- x$file == i
    lines(x$rt[ii], x$intensity[ii], col = col[diffPheno[match(i, unifile)]])
  }
  
  legend("topright", col = col, legend = names(col), lty = 1, bty = "n", lwd = 3)
}