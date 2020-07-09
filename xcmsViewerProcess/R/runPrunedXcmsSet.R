#' mapping ms2 chromPeaksMS2
#' @param x the MSnExp oject
#' @param mtab_files paths of meta table, often returned by function \code{peakPicking }
#' @param itab_files paths of the intensity table, often returned by function \code{peakPicking }
#' @param fun_parallel the parallel function, could be \code{mclapply} or \code{bplapply}
#' @param ... the parameters passed to parallel function
#' @import parallel
#' @import BiocParallel
#' @noRd
chromPeaksMS2 <- function(x, mtab_files, itab_files, fun_parallel = parallel::mclapply, ...) {
  
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
  cpeaks
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
#'  The elements in the returned list could be passed to \code{defineFeatures}
#' @export
#' 
peakPicking <- function(
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
  
  flt <- postfilterParam
  flt$postfilter <- NULL
  flt$BPPARAM <- NULL
  
  list(itab = file_itab, 
       mtab = file_mtab, 
       XCMSnExp = file_XCMSnExp, 
       xcmsScanFilter = flt, 
       peakPickingParam = param, 
       postFilter = postfilterParam$postfilter)
}

#' define features for individual pipeline 
#' @param files the file paths of XCMSnExp
#' @param mtab_files the file path of meta tables, often returned by function \code{peakPicking}
#' @param pgParam the peak grouping parameter, passed to \code{groupChromPeaks}
#' @param rtParam the retention time adjustment parameters, passed to \code{adjustRtime}
#' @export
#' @importFrom stringr str_split_fixed
#' @return It returns an object of class \code{XCMSnExp}
#' 
defineFeatures <- function(files, mtab_files, rtParam, pgParam) {
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


#' xcms summarize function for the individual pipeline
#' @param x the MSnExp oject
#' @param mtab_files paths of meta table, often returned by function \code{peakPicking }
#' @param itab_files paths of the intensity table, often returned by function \code{peakPicking }
#' @param fun_parallel the parallel function, could be \code{mclapply} or \code{bplapply}
#' @param ... the parameters passed to parallel function
#' @import parallel
#' @import BiocParallel
#' @export
#' @return An object of class "prunedXcmsSet"
#' 
runPrunedXcmsSet <- function( 
  
  # peakPicking
  files, 
  peakPickingParam = MatchedFilterParam(fwhm = 7.5), 
  tmpdir="xcmsViewerTemp",  
  postfilterParam = list(
    postfilter = c(5, 1000), 
    ms1.noise = 100, ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 20,
    ms2.noise = 30, ms2.maxPeaks = 100, ms2.maxIdenticalInt = 6, 
    BPPARAM=bpparam()),
  
  # phenotype data
  pheno = NULL,
  
  # feature identification
  retentionTimeParam = PeakGroupsParam(
    minFraction = 0.5,
    extraPeaks = 1,
    smooth = "loess",
    span = 0.5,
    family = "gaussian",
    peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
    subset = integer(),
    subsetAdjust = c("average", "previous")), 
  PeakGroupsParam = PeakDensityParam(rep(1, 4)),
  fun_parallel = parallel::mclapply, 
  ... 
  ) {
  
  # peakPicking and filter
  pp <- peakPicking(
    files = files, param = peakPickingParam, tmpdir=tmpdir, postfilterParam = postfilterParam
    )
  
  # feature identification, the RT will also be identified
  df <- defineFeatures(
    files = pp$XCMSnExp, mtab_files = pp$mtab, rtParam = retentionTimeParam,  pgParam = PeakGroupsParam
    )
  
  cat("Extracting extended chrom peaks ...\n")
  peaks <- chromPeaksMS2(
    df, mtab_files = pp$mtab, itab_files = pp$itab, fun_parallel = fun_parallel, ...
  )
  peaks$ID <- rownames(peaks)
  cat(" done! \n")
  
  cat("Extracting extended feature table ...")
  fd <- as.data.frame(featureDefinitions(df))
  fd$ID <- rownames(fd)
  fd$annot_ms1 <- character(nrow(fd))
  fd$annot_ms2 <- character(nrow(fd))
  fd <- fd[, .xcmsViewerInternalObjects$xcmsFeatureSet_fdata_name]
  peaks$sample <- as.factor(peaks$sample)
  fv <- sapply(fd$peakidx, function(i) {
    qm <- peaks[i, ]
    qm <- qm[order(qm$into, decreasing = TRUE), ]
    list(
      intensity = tapply(qm$into, qm$sample, "[", 1),
      feat = tapply(rownames(qm), qm$sample, "[", 1)
    )
  })
  peaks$sample <- as.integer(as.character(peaks$sample))
  
  cat("Summarizing experiment ...")
  if (is.null(pheno)) {
    pheno <- pData(df)
    pheno$num <- 1:nrow(pheno)
  }
  
  int <- do.call(rbind, fv[1, ])
  masterPeak <- do.call(rbind, fv[2, ])
  rownames(int) <- rownames(masterPeak) <- rownames(fd)
  colnames(int) <- colnames(masterPeak) <- rownames(pheno)
  
  peaks$masterPeak <- ""
  peaks$masterPeak[fastmatch::"%fin%"(peaks$ID, masterPeak)] <- "+"
  
  new("xcmsPeak", 
      table = peaks[, .xcmsViewerInternalObjects$xcmsPeak_table_column],
      param = pp$peakPickingParam,
      postFilter = pp$postFilter)
  
  
  obj_xcmsPeak <- obj_xcmsFeatureSet <- new(
    "xcmsFeatureSet", 
    intensity = int, masterPeak = masterPeak, 
    phenoData = AnnotatedDataFrame(pheno),
    featureData = AnnotatedDataFrame(fd))
  
  itable <- do.call(rbind, lapply(pp$itab, readRDS))
  mtable <- do.call(rbind, lapply(pp$mtab, readRDS))
  obj_xcmsScan <- new("xcmsScan", meta = mtable, intensity = itable, filter = pp$xcmsScanFilter)

  res <- list(
    xcmsFeatureSet = obj_xcmsFeatureSet,
    xcmsPeak = obj_xcmsPeak,
    xcmsScan = obj_xcmsScan
  )
  res
}

