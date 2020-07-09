
#' xcms annotation function for the individual pipeline
#' @param x an object of XCMSnExp
#' @param mode mode. pos or neg
#' @param massTab a data.frame of mass table having at least one column named as "monoisotopic_molecular_weight"
#' @param refSpectra reference MS2 spectra used annotate the experimental MS2 spectra
#' @param fun_parallel the parallel function, could be \code{mclapply} or \code{bplapply}
#' @param ... the parameters passed to parallel function
#' @import parallel
#' @import BiocParallel
#' @export
#' 
xcmsAnnotate <- function(
  x, mode = c("pos", "neg")[1], massTab=NULL, refSpectra = NULL, fun_parallel = mclapply, ...
) {
  tmp <- list(...)
  features <- x$features
  features$meta$Annotation <- ""
  peaks <- x$peaks
  pheno <- x$pheno
  
  itable <- do.call(rbind, lapply(x$scanIntensityTab, readRDS))
  mtable <- do.call(rbind, lapply(x$scanMetaTab, readRDS))
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
    try(.validate_statsXCMS(res))
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
    adduct = as.character(at$name),
    massdiff = at$massdiff,
    nmol = as.numeric(at$nmol),
    ips = as.numeric(at$ips),
    stringsAsFactors = FALSE)
  
  ll <- fun_parallel(1:nrow(features$meta), function(i) {
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
  res$annotationMass <- ll
  
  res$features$meta$Annotation <- sapply(ll, function(x) {      
    i <- 1:min(5, nrow(x))
    paste(x$name[i], collapse = ";")
  })
  
  if (is.null(refSpectra)) {
    try(.validate_statsXCMS(res))
    return(res)
  }
  
  cat("Annotating MS2 spectra ... \n")
  ms2a <- .ms2map(itab_files = x$scanIntensityTab, 
                  featureMeta = features$meta, 
                  valideRefs = ll, 
                  mode = mode, 
                  refs = refSpectra, 
                  peaks = peaks,
                  fun_parallel = fun_parallel,
                  ...)
  
  if (is.null(ms2a)) {
    try(.validate_statsXCMS(res))
    return(res)
  }
  
  ap <- unique(ms2a$annot_peaks)
  ail <- list(
    meta = refSpectra$meta[refSpectra$meta$id2 %fin% ap, ],
    peakList = refSpectra$peakList[ap]
  )
  colnames(ail$meta) <- gsub("-", "_", colnames(ail$meta))
  
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
    # s_massAt <- cbind(s_massAt, .dd)
    s_massAt$annot_peaks <- .dd$annot_peaks
    s_massAt$query_peaks <- .dd$query_peaks
    s_massAt$cos <- .dd$cos
    s_massAt$database_id <- .dd$database_id
    s_massAt$atScore <- s_massAt$score
    .cv <- !is.na(s_massAt$cos)
    s_massAt$atScore[.cv] <- s_massAt$atScore[.cv] + s_massAt$cos[.cv]
    s_massAt[order(s_massAt$atScore, decreasing = TRUE), ]
  })
  names(vx) <- features$meta$ID
  ll <- vx
  
  features$meta$Annotation <- sapply(vx, function(x) {      
    i <- 1:min(5, nrow(x))
    i <- paste(x$name[i], collapse = ";")
    if (nchar(i) > 300)
      i <- paste0(substr(i, 1, 300), "...")
    i
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
  try(.validate_statsXCMS(res))
  res
}

#' Individual pipeline for summarizeExp
#' @description the individual pipeline process each file separately so less memory
#'   used, but slower.
#' @param x a data.frame where the first columns should be named as "file". It 
#'   contains the full path of mzXML or mzML file
#' @param mode mode. pos or neg
#' @param peakPickingParam peak picking parameters, passed to \code{findChromPeaks}
#' @param postfilterParam passed to \code{chromPeaksPostFilter}
#' @param PeakGroupsParam the peak grouping parameter, passed to \code{groupChromPeaks}
#' @param retentionTimeParam the retention time adjustment parameters, passed to \code{adjustRtime}
#' @param massTab a data.frame of mass table having at least one column named as "monoisotopic_molecular_weight"
#' @param refSpectra reference MS2 spectra used annotate the experimental MS2 spectra
#' @param tmpdir the temporary directory to store the processed results
#' @param parallelFun the parallel function, could be \code{mclapply} or \code{bplapply}
#' @param ... the parameters passed to parallel function
#' @export
#' 
summarizeExp <- function(
  x,
  mode = c("pos", "neg")[1],
  peakPickingParam = MatchedFilterParam(fwhm = 7.5),
  postfilterParam = list(
    postfilter = c(5, 1000), 
    ms1.noise = 100, ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 20,
    ms2.noise = 30, ms2.maxPeaks = 100, ms2.maxIdenticalInt = 6, 
    BPPARAM=bpparam()),
  PeakGroupsParam = PeakGroupsParam(
    minFraction = 0.5,
    extraPeaks = 1,
    smooth = "loess",
    span = 0.5,
    family = "gaussian",
    peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
    subset = integer(),
    subsetAdjust = c("average", "previous")),
  retentionTimeParam = NearestPeaksParam(),
  massTab = NULL,
  refSpectra = NULL, 
  parallelFun = bplapply, 
  tmpdir = NA,
  ...
) {
  
  if (is.na(tmpdir)) 
    tmpdir <- "xcmsViewerTemp"
  
  files <- x$file
  x$file <- basename(x$file)
  
  print(Sys.time())
  cat("Peak picking ...\n")
  ff <- peakPicking(
    files = files,
    param = peakPickingParam,
    postfilterParam = postfilterParam,
    tmpdir = tmpdir
  )
  
  print(Sys.time())
  cat("Feature defining ...\n")
  xd <- defineFeatures(
    files = ff$XCMSnExp,
    mtab_files = ff$mtab,
    rtParam = retentionTimeParam,
    pgParam = PeakGroupsParam )
  saveRDS(xd, file = file.path(tmpdir, "04_defineFeatures.RDS"))
  
  print(Sys.time())
  cat("Creating xcmsSummarize object ...\n")
  xdata <- xcmsSummarize(
    xd, mtab_files = ff$mtab, itab_files = ff$itab, fun_parallel = parallelFun, ...
  )
  xdata$pheno <- x
  saveRDS(xdata, file = file.path(tmpdir, "05_xcmsSummarize.RDS"))
  
  print(Sys.time())
  cat("Annotation ...\n")
  st <- xcmsAnnotate(
    x = xdata, massTab = massTab, refSpectra = refSpectra, mode = mode, 
    fun_parallel = parallelFun, ...
  )
  saveRDS(st, file = file.path(tmpdir, "xcmsViewerData.RDS"))
  
}

#' Internal function for ms2mapping
#' @param itab_files file paths of intensity table
#' @param featureMeta feature meta data
#' @param mode pos or neg
#' @param valideRefs valideRefs
#' @param refs refs
#' @param peaks the peaks
#' @param fun_parallel the parallel function
#' @param ... other parameters passed to parallel
#' @noRd
#' 
.ms2map <- function(
  itab_files, featureMeta, valideRefs, mode, refs, peaks, fun_parallel, ...
) {
  
  ms2an <- list()
  featureId <- featureMeta$ID
  featurePeakIds <- featureMeta$peakidx
  names(featurePeakIds) <- featureId
  
  for (i in seq_along(itab_files)) {
    fn <- itab_files[i]
    itable <- readRDS(fn)
    print(sprintf("%s/%s", i, length(itab_files)))
    ms2an[[i]] <- fun_parallel( featureId, function(j) {
      spectraCosDist(
        ftid = j, peakidx = featurePeakIds[[j]],
        possibleRefs = valideRefs, mode = mode, refSpectra = refs,
        peakTab = peaks, itab = itable
      )
    }, ...)
  }
  ms2an <- unlist(ms2an, recursive = FALSE)
  ms2a <- do.call(rbind, ms2an)
  if (is.null(ms2a))
    return(NULL)
  colnames(ms2a) <- gsub("-", "_", colnames(ms2a))
  ms2a
}
