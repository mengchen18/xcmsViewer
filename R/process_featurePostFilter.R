#' post filter of features
#' @param x An XCMSnExp object after peak calling
#' @param nComplete the number of complete samples in a 
#'  group (is length = 1) /all groups (if length = group numbers). If can be given
#'  integers (>= 1) or the percent of groups size (0, 1).#' 
#' @export

featurePostfilter <- function(x, nComplete = 1) {
  fd <- featureDefinitions(x)
  gpc <- setdiff(colnames(fd), c("mzmed", "mzmin", "mzmax", "rtmed", 
                                 "rtmin", "rtmax", "npeaks", "peakidx"))
  
  if (!is.numeric(nComplete))
    stop("nComplete should numeric!")
  
  if (length(nComplete) == 1) {
    
    if (nComplete >= 1 || nComplete == 0) {
      ct <- structure(rep(nComplete, length(gpc)), names = gpc)
    } else if (nComplete < 1 & nComplete > 0) {
      ct <- structure(sapply(fd[gpc], max)*nComplete, names = gpc)
    } else
      stop("nComplete should >= 0")

    ll <- lapply(gpc, function(i) {
      fd[[i]]  >= ct[i]
    })
    f <- Reduce("|", ll)
    
  } else if (length(nComplete) == length(gpc)) {
    
    if (!all(names(nComplete) %in% gpc))
      stop("If nComplete has the same length as group number, it should be a named
           vector. The names are the same as group names.")
    
    if (all(nComplete >= 1 | nComplete == 0)) {
      ct <- nComplete[gpc]
    } else if (all(nComplete < 1 & nComplete > 0)) {
      ct <- sapply(fd[gpc], max)*nComplete[gpc]
    } else
      stop("nComplete should all >= 1 or within (0, 1)!")

    ll <- lapply(gpc, function(i) {
      fd[[i]]  >= ct[i]
    })
    f <- Reduce("&", ll)
    
  } else 
    stop("length of nComplete should be either 1 or the same as group number!")
  

  cat(sprintf("Removing %s features ...\n", sum(!f)))
  peaks <- chromPeaks(x)
  fdf <- fd[f, ]
  idl <- lapply(fdf$peakidx, function(x) rownames(peaks)[x])
  ii <- rownames(peaks) %in% unlist(idl)
  peaks <- peaks[ii, ]
  fdf$peakidx <- lapply(idl, fastmatch::fmatch, rownames(peaks))
  rownames(peaks) <- xcms:::.featureIDs(nrow(peaks), prefix = "CP")
  
  chromPeaks(x) <- peaks
  featureDefinitions(x) <- fdf
  x
}

