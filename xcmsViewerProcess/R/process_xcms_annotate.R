#' @describeIn summarizeExp Annotating the summarized xcms object
#' @export

xcms_annotate <- function(x, mode = c("pos", "neg")[1], massTab=NULL, refSpectra = NULL, ...) {
  
  features <- x$features
  features$meta$Annotation <- ""
  peaks <- x$peaks
  mtable <- x$scanMetaTab
  itable <- x$scanIntensityTab
  pheno <- x$pheno
  
  res <- list(
    features = features,
    peaks = peaks,
    scanMetaTab = mtable,
    scanIntensityTab = itable,
    annotationMass = NULL,
    annotationFragment = NULL,
    matchedRefFragments = NULL,
    pheno = pheno
  )
  
  if (is.null(massTab)) {
    .validate_statsXCMS(res)
    return(res)
  }
  
  cat("Annotating using mass ... \n")
  maitEnv <- environment()
  data("MAITtables", package = "MAIT", envir = maitEnv)
  if (mode == "pos") {
    at = maitEnv$posAdducts
  } else if (mode == "neg") 
    at = maitEnv$negAdducts else 
      stop("'mode' should be either pos or neg!")
  at <- data.frame(
    adduct = as.character(at$name)[at$nmol==1],
    massdiff = at$massdiff[at$nmol==1],
    stringsAsFactors = FALSE)
  
  ll <- bplapply(1:nrow(features$meta), function(i) {
    r <- massQuery(m = features$meta$mzmed[i], tolppm = 10, refTab = massTab, 
              addTable = at, ID = rownames(features$meta)[i]) 
    if (!is.null(r)) {
      r$annot_peaks <- as.character(NA)
      r$query_peaks <- as.character(NA) 
      r$cos <- as.numeric(NA)
      r$database_id <- as.character(NA)
      r$atScore <- as.numeric(NA)
      }    
    r
  }, ...)
  names(ll) <- features$meta$ID
  gc()
  cat("done!\n")
  res$annotationMass <- ll

  res$features$meta$Annotation <- sapply(ll, function(x) {      
    i <- 1:min(5, nrow(x))
    paste(x$name[i], collapse = ";")
  })
  
  if (is.null(refSpectra)) {
    .validate_statsXCMS(res)
    return(res)
  }
  
  cat("Annotating MS2 spectra ... \n")
  ms2an <- bplapply( features$meta$ID, function(x) {
    spectraCosDist(
      ftid = x, peakidx = features$meta[x, "peakidx"][[1]], 
      possibleRefs = ll, mode = mode, refSpectra = refSpectra, 
      peakTab = peaks, itab = itable
    )
  }, ...)
  names(ms2an) <- features$meta$ID
  
  if (all(sapply(ms2an, is.null))) {
    .validate_statsXCMS(res)
    return(res)
  }
  
  ms2a <- do.call(rbind, ms2an)
  colnames(ms2a) <- gsub("-", "_", colnames(ms2a))
  ap <- unique(ms2a$annot_peaks)
  ail <- list(
    meta = refSpectra$meta[refSpectra$meta$id2 %fin% ap, ],
    peakList = refSpectra$peakList[ap]
  )
  colnames(ail$meta) <- gsub("-", "_", colnames(ail$meta))
  gc()
  cat("done!\n")
  
  cat("Combining MS2 and MS1 annotations ... \n")
  vx <- lapply(features$meta$ID, function(id) {
    ## annotation by mass comparison
    s_massAt <- ll[[id]]
    if (is.null(s_massAt))
      return(NULL)
    i <- fmatch(id, features$meta$ID)
    s_peakids <- features$meta[i, "peakidx"][[1]] # get peak ids
    s_ms2scans <- peaks[features$meta[i, "peakidx"][[1]], "ms2Scan"] # get MS2 scans
    s_ms2scans <- unlist(strsplit(s_ms2scans, split = ";"))
    s_ms2scans_meta <- mtable[ mtable$ID %fin% s_ms2scans, ]
    s_ms2match <- ms2a[ms2a$query_peaks %fin% s_ms2scans_meta$ID, ]
    s_ms2match <- s_ms2match[order(s_ms2match$cos, decreasing = TRUE), ]
    s_ms2match <- s_ms2match[!duplicated(s_ms2match$database_id), ]
    
    .dd <- data.frame(
      annot_peaks = rep(NA, nrow(s_massAt)),
      query_peaks = NA,
      cos = NA,
      database_id = NA)
    
    for (ii in 1:nrow(s_ms2match)) {
      i <- which(s_ms2match$database_id[ii] == s_massAt$accession)
      .dd[i, ] <- s_ms2match[ii, ]
    }
    if (is.logical(.dd$cos)) {
      class(.dd$annot_peaks) <- "character"
      class(.dd$query_peaks) <- "character"
      class(.dd$cos) <- "numeric"
      class(.dd$database_id) <- "character"
    }
    s_massAt <- cbind(s_massAt, .dd)
    s_massAt$atScore <- s_massAt$score
    .cv <- !is.na(s_massAt$cos)
    s_massAt$atScore[.cv] <- s_massAt$atScore[.cv] + s_massAt$cos[.cv]
    s_massAt[order(s_massAt$atScore, decreasing = TRUE), ]
  })
  names(vx) <- features$meta$ID
  ll <- vx
  
  features$meta$Annotation <- sapply(vx, function(x) {      
    i <- 1:min(5, nrow(x))
    paste(x$name[i], collapse = ";")
  })
  
  res <- list(
    features = features,
    peaks = peaks,
    scanMetaTab = mtable,
    scanIntensityTab = itable,
    annotationMass = ll,
    annotationFragment = ms2a,
    matchedRefFragments = ail,
    pheno = pheno
  )
  .validate_statsXCMS(res)
  res
}