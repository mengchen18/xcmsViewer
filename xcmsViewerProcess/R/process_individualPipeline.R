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
#' @param fun_parallel the parallel function, could be \code{mclapply} or \code{bplapply}
#' @param ... the parameters passed to parallel function
#' @export
#' 
indiv_summarizeExp <- function(
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
  ff <- indiv_peakPicking(
    files = files,
    param = peakPickingParam,
    postfilterParam = postfilterParam,
    tmpdir = tmpdir,
    BPPARAM = BPPARAM
  )
  
  print(Sys.time())
  cat("Feature defining ...\n")
  xd <- indiv_defineFeatures(
    files = ff$XCMSnExp,
    mtab_files = ff$mtab,
    rtParam = retentionTimeParam,
    pgParam = PeakGroupsParam )
  saveRDS(xd, file = file.path(tmpdir, "04_defineFeatures.RDS"))
  
  print(Sys.time())
  cat("Creating xcmsSummarize object ...\n")
  xdata <- indiv_xcmsSummarize(
    xd, mtab_files = ff$mtab, itab_files = ff$itab, fun_parallel = parallelFun, ...
  )
  data$pheno <- x
  saveRDS(xdata, file = file.path(tmpdir, "05_xcmsSummarize.RDS"))
  
  print(Sys.time())
  cat("Annotation ...\n")
  st <- indiv_xcmsAnnotate(
    x = xdata, massTab = massTab, refSpectra = refSpectra, mode = mode, 
    fun_parallel = parallelFun, ...
  )
  
}





#' peak picking for individual pipeline 
#' @description peaks picking are peroformed one raw file after another so that less memory
#'   used. 
#' @param files a character vector stores the the mzXML or mzML files
#' @param param peak picking parameters, passed to \code{findChromPeaks}
#' @param tmpdir the temporary directory to store the processed results
#' @param postfilterParam passed to \code{chromPeaksPostFilter}
#' @return it returns a list of three elements:
#'  1. itab = file paths of intensity tables
#'  2. mtab = file  paths of the meta tables
#'  3. XCMSnExp = file paths of XCMSnExp
#'  The elements in the returned list could be passed to \code{indiv_defineFeatures}
#' @export
#' 
indiv_peakPicking <- function(
  files, param, tmpdir="xcmsViewerTemp",  
  postfilterParam = list(
    postfilter = c(5, 1000), 
    ms1.noise = 100, ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 20,
    ms2.noise = 30, ms2.maxPeaks = 100, ms2.maxIdenticalInt = 6, 
    BPPARAM=bpparam())
) {
  
  # creating temp dirs 
  if (!file.exists(tmpdir))
    dir.create(tmpdir) else {
      unlink(tmpdir, recursive = TRUE)
      dir.create(tmpdir)
    }
  
  dir_itab <- file.path(tmpdir, "01_itab")
  dir_mtab <- file.path(tmpdir, "02_mtab")
  dir_peak <- file.path(tmpdir, "03_XCMSnExp_peakPicking")
  
  dir.create(dir_itab)
  dir.create(dir_mtab)
  dir.create(dir_peak)
  
  file_itab <- c()
  file_mtab <- c()
  file_XCMSnExp <- c()
  
  # loop
  loop <- 0
  nfiles <- length(files)
  nd <- nchar(nfiles)
  for (f1 in files) {
    loop <- loop + 1
    print(sprintf("%s/%s", loop, nfiles))
    
    fid <- paste0("F", formatC(loop, width=nd, flag="0"), ".")
    f1b <- gsub("\\..*$", ".Rds", basename(f1))
    
    d1 <- MSnbase::readMSData(f1, mode = "onDisk")
    xdata <- xcms::findChromPeaks( d1, param = param, msLevel = 1L)
    postfilterParam$x <- xdata
    xdata <- do.call(chromPeaksPostFilter, args = postfilterParam)
    
    itab <- xdata$itable
    itab$ID <- sub("^F1.", fid, itab$ID)
    saveRDS(itab, file = file.path(dir_itab, f1b))
    file_itab <- c( file_itab, file.path(dir_itab, f1b) )
    
    mtab <- xdata$mtable
    mtab$fromFile <- loop
    mtab$ID <- sub("F1.", fid, mtab$ID)
    saveRDS(mtab, file = file.path(dir_mtab, f1b))
    file_mtab <- c( file_mtab, file.path(dir_mtab, f1b) )
    
    saveRDS(xdata$XCMSnExp, file = file.path(dir_peak, f1b))
    file_XCMSnExp <- c( file_XCMSnExp, file.path(dir_peak, f1b) )
  }
  list(itab = file_itab, mtab = file_mtab, XCMSnExp = file_XCMSnExp)
}

#' define features for individual pipeline 
#' @param files the file paths of XCMSnExp
#' @param mtab_files the file path of meta tables, often returned by function \code{indiv_peakPicking}
#' @param pgParam the peak grouping parameter, passed to \code{groupChromPeaks}
#' @param rtParam the retention time adjustment parameters, passed to \code{adjustRtime}
#' @export
#' @return It returns an object of class \code{XCMSnExp}
#' 
indiv_defineFeatures <- function(files, mtab_files, rtParam, pgParam) {
  exps <- lapply(files, readRDS)
  xdata <- do.call(c, exps)
  xdata <- groupChromPeaks(xdata, param = pgParam)
  xdata <- adjustRtime(xdata, param = rtParam)
  object <- groupChromPeaks(xdata, param = pgParam)
  
  if (hasAdjustedRtime( object )) {
    rt <- rtime(object)
    id <- stringr::str_split_fixed(names(rt),  pattern = "\\.", n = 2)
    idu <- unique(id[, 1])
    
    for (i in 1:length(mtab_files)) {
      x <- readRDS(mtab_files[i])
      ii <- id[, 1] == idu[i]
      x$rt.unadjusted <- x$rt
      x$rt <- rt[ii]
      saveRDS(x, file = mtab_files[i])
    }
  }
  
  object
}


#' mapping ms2 chromPeaksMS2
#' @param x the MSnExp oject
#' @param mtab_files paths of meta table, often returned by function \code{indiv_peakPicking }
#' @param itab_files paths of the intensity table, often returned by function \code{indiv_peakPicking }
#' @param fun_parallel the parallel function, could be \code{mclapply} or \code{bplapply}
#' @param ... the parameters passed to parallel function
#' @import parallel
#' @import BiocParallel
#' 
.indiv_chromPeaksMS2 <- function(x, mtab_files, itab_files, fun_parallel = parallel::mclapply, ...) {
  
  if (!inherits(x, "MSnExp"))
    stop("x should be an object of class MSnExp!")
  
  cpeaks <- xcms::chromPeaks(x)
  cpeaks <- as.data.frame(cpeaks)
  
  cat("Mapping MS2 spectra ...\n")
  t0 <- max(cpeaks$sample)
  l <- lapply(unique(cpeaks$sample), function(fromFile) {
    print(sprintf("%s/%s", fromFile, t0))
    cpsub <- cpeaks[cpeaks$sample == fromFile, ]
    mtab <- readRDS(mtab_files[fromFile])
    mtab <- mtab[mtab$msLevel > 1, ]
    fun_parallel(1:nrow(cpsub), function(i) {
      c1 <- cpsub[i, ]
      ic <- which(
        dplyr::between(mtab[["rt"]], c1[["rtmin"]], c1[["rtmax"]] ) & 
          dplyr::between(mtab$precMz, c1[["mzmin"]], c1[["mzmax"]])
      )
      paste(mtab$ID[ic], collapse = ";")
    }, ...)
  })
  cpeaks$ms2Scan <- unlist(l)
  cpeaks$QC <- "+"
  cpeaks
}

#' xcms summarize function for the individual pipeline
#' @param x the MSnExp oject
#' @param mtab_files paths of meta table, often returned by function \code{indiv_peakPicking }
#' @param itab_files paths of the intensity table, often returned by function \code{indiv_peakPicking }
#' @param fun_parallel the parallel function, could be \code{mclapply} or \code{bplapply}
#' @param ... the parameters passed to parallel function
#' @import parallel
#' @import BiocParallel
indiv_xcmsSummarize <- function( x, mtab_files, itab_files, fun_parallel = mclapply, ... ) {
  
  cat("Extracting extended chrom peaks ...\n")
  peaks <- .indiv_ChromPeaksMS2(
    x, mtab_files = mtab_files, itab_files = itab_files, fun_parallel = fun_parallel, ...
  )
  
  peaks$sample <- as.factor(peaks$sample)
  peaks$ID <- rownames(peaks)
  cat(" done! \n")
  
  pheno <- pData(x)
  pheno$num <- 1:nrow(pheno)
  
  cat("Extracting extended feature table ...")
  fd <- as.data.frame(featureDefinitions(x))
  fd$ID <- rownames(fd)
  fv <- sapply(fd$peakidx, function(i) {
    qm <- peaks[i, ]
    qm <- qm[!is.na(qm$QC), ]
    qm <- qm[order(qm$into, decreasing = TRUE), ]
    list(
      intensity = tapply(qm$into, qm$sample, "[", 1),
      feat = tapply(rownames(qm), qm$sample, "[", 1)
    )
  })
  
  features <- list(
    meta = fd,
    intensities = data.frame(
      ID = rownames(fd), 
      do.call(rbind, fv[1, ]),
      stringsAsFactors = FALSE),
    masterPeaks = data.frame(
      ID = rownames(fd), 
      do.call(rbind, fv[2, ]),
      stringsAsFactors = FALSE)
  )
  rownames(features$intensities) <- rownames(features$masterPeaks) <- rownames(features$meta)
  features$meta$QC <- apply(features$masterPeaks, 1, function(x) {
    ids <- unlist(strsplit(x, ";"))
    r <- ""
    if (any(peaks$QC[peaks$ID %fin% ids] == "+"))
      r <- "+"
    r
  })
  
  v <- lapply(features$masterPeaks[setdiff(colnames(features$masterPeaks), "ID")], function(x)
    na.omit(setdiff(unique(unlist(strsplit(x, ";"))), ""))
  )
  vv <- unlist(v)
  peaks$masterPeak <- ""
  peaks$masterPeak[peaks$ID %fin% vv] <- "+"
  cc <- c("ID",  "QC", "sample", "masterPeak", "into", 
          "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "intb", "intf", "maxo", "sn", 
          "ms2Scan", "rsq", "rtgap", "intgap", "rtintgap", "truncated", "b")
  cc <- cc[cc %in% colnames(peaks)]
  peaks <- peaks[, cc]
  cat(" done! \n")
  
  res <- list(
    features = features,
    peaks = peaks,
    scanMetaTab = mtab_files,
    scanIntensityTab = itab_files,
    pheno = pheno
  )
  res
}


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
indiv_xcmsAnnotate <- function(
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
  ms2a <- .indiv_ms2map(itab_files = x$scanIntensityTab, 
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

#' Internal function for ms2mapping
.indiv_ms2map <- function(
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
  # # names(ms2an) <- features$meta$ID
  ms2an <- unlist(ms2an, recursive = FALSE)
  ms2a <- do.call(rbind, ms2an)
  if (is.null(ms2a))
    return(NULL)
  colnames(ms2a) <- gsub("-", "_", colnames(ms2a))
  ms2a
}

