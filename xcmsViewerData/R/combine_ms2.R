#' combine ms2 annotation object
#' @param ... the annotation object. Each should be a list object containing at least two elements named 'meta' and 'peakList'
#' @param list a list of annotation object
#' @export
#' 
combine_ms2 <- function(..., list = NULL) {
  l <- list(...)
  l <- c(l, list)
  for (i in l) {
    if (!all(names(i) %in% c("meta", "peakList")))
      stop("All ms2 list should contain at least two elements named 'meta' and 'peakList'!")
  }
  l_meta <- lapply(l, "[[", "meta")
  l_peakList <- lapply(l, "[[", "peakList")
  
  list(
    meta = do.call(rbind, l_meta),
    peakList = unlist(l_peakList, recursive = FALSE)
  )
}

