#' function to extract ms1 XICs given a simple intensity table format
#' @param itab an intensity table
#' @param mtab the spectra meta table of the experiment
#' @param rt range of retention time
#' @param mz range of mz
#' @param file which file to extract
#' @param na.rm logical, whether remove NA values
#' @importFrom dplyr between
#' @importFrom fastmatch %fin%
#' @importFrom fastmatch fmatch
#' @importFrom stats na.omit
#' @import graphics
#' @export
#' @examples 
#' x <- readRDS(system.file(package="xcmsViewerApp", "extdata", "00_exampleData_processed.RDS"))
#' mtable <- x$scanMetaTab
#' itable <- x$scanIntensityTab
#' chrom <- eic(itab=itable, mtab=mtable, rt = c(0, Inf), mz=c(225, 255), file=1, na.rm = TRUE)

eic <- function(itab, mtab, rt, mz, file, na.rm = TRUE) {
  
  if (missing(file))
    file <- unique(mtab$fromFile)
  
  irt <- which(dplyr::between(mtab$rt, rt[1], rt[2]) & mtab$msLevel == 1 & mtab$fromFile %fin% file)
  mtab_subset <- mtab[irt, ]
  imz <- which(itab$ID %fin% mtab_subset$ID & dplyr::between(itab$mz, mz[1], mz[2]))
  if (length(imz) == 0)
    return()
  
  l <- itab[imz, ]
  intensity <- tapply(l$intensity, l$ID, sum)
  it <- fmatch(names(intensity), mtab$ID)
  l <- data.frame(
    rt = mtab$rt[it],
    intensity = intensity,
    file = mtab$fromFile[it]
  )
  if (na.rm)
    l <- na.omit(l)
  if (length(missingfile <- setdiff(file, l$file)) > 0)
    l <- rbind(l,
      data.frame(
        rt = min(l$rt),
        intensity = min(l$intensity),
        file = missingfile
        )
      )
    l <- l[order(l$file), ]
  l
}
