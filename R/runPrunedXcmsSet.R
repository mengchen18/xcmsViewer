#' mapping ms2 chromPeaksMS2
#' @param x the MSnExp object
#' @param mtab_files paths of meta table, often returned by function \code{peakPicking }
#' @param fun_parallel the parallel function, could be \code{mclapply} or \code{bplapply}
#' @param ... the parameters passed to parallel function
#' @import parallel
#' @import BiocParallel
#' @noRd
#' 
chromPeaksMS2 <- function(x, mtab_files, fun_parallel = parallel::mclapply, ...) {
  
  if (!inherits(x, "MSnExp"))
    stop("x should be an object of class MSnExp!")
  
  cpeaks <- xcms::chromPeaks(x)  
  cat("Mapping MS2 spectra ...\n")
  t0 <- max(cpeaks[, "sample"])
  l <- lapply(unique(cpeaks[, "sample"]), function(fromFile) {
    cpsub <- cpeaks[cpeaks[, "sample"] == fromFile, ]
    mtab <- readRDS(mtab_files[fromFile])
    mtab <- mtab[mtab$msLevel > 1, ]
    fun_parallel(1:nrow(cpsub), function(i) {
      c1 <- cpsub[i, ]
      ic <- which(
        dplyr::between(mtab[["rt"]], c1[["rtmin"]], c1[["rtmax"]] ) & 
          dplyr::between(mtab$precMz, c1[["mzmin"]], c1[["mzmax"]])
      )
      list(paste(mtab$ID[ic], collapse = ";"), any(mtab$validMS2[ic]))
    }, ...)
  })
  v <- character(nrow(cpeaks))
  vms2 <- rep(FALSE, nrow(cpeaks))
  for (fromFile in unique(cpeaks[, "sample"])) {
    i0 <- cpeaks[, "sample"] == fromFile
    v[i0] <- sapply(l[[fromFile]], "[[", 1)
    vms2[i0] <- sapply(l[[fromFile]], "[[", 2)
  }
  data.frame(
    ms2Scan = v,
    validMS2 = vms2
    )
}


#' peak picking for individual pipeline 
#' @description peaks picking are peroformed one raw file after another so that less memory
#'   used. 
#' @param files a character vector stores the the mzXML or mzML files
#' @param param peak picking parameters, passed to \code{findChromPeaks}
#' @param tmpdir the temporary directory to store the processed results
#' @param postfilter passed to \code{chromPeaksPostFilter}
#' @param noisefilter a list, passed to \code{getMetaIntensityTable}
#' @param BPPARAM arguments passed to chromPeaksPostFilter, which further pass the argument to bplapply
#' @return it returns a list of three elements:
#'  1. itab = file paths of intensity tables
#'  2. mtab = file  paths of the meta tables
#'  3. XCMSnExp = file paths of XCMSnExp
#'  The elements in the returned list could be passed to \code{defineFeatures}
#' @export
#' 
peakPicking <- function(
  files, param, tmpdir="xcmsViewerTemp",  
  postfilter = c(5, 1000),
  noisefilter = list(
    ms1.noise = 100, ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 20,
    ms2.noise = 30, ms2.maxPeaks = 100, ms2.maxIdenticalInt = 6),
  BPPARAM=bpparam()
) {
  
  # creating temp dirs 
  if (!file.exists(tmpdir))
    dir.create(tmpdir) else {
      if (length(list.files(tmpdir)) > 0)
        stop("tmpdir is not empty, this may lead to data loss! Please select an empty folder!")
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
    noisefilter$x <- xdata
    mitab <- do.call(getMetaIntensityTable, noisefilter)
    
    itab <- mitab$itab
    itab$ID <- sub("^F1.", fid, itab$ID)
    saveRDS(itab, file = file.path(dir_itab, f1b))
    file_itab <- c( file_itab, file.path(dir_itab, f1b) )
    
    mtab <- mitab$mtab
    mtab$fromFile <- loop
    mtab$ID <- sub("F1.", fid, mtab$ID)
    saveRDS(mtab, file = file.path(dir_mtab, f1b))
    file_mtab <- c( file_mtab, file.path(dir_mtab, f1b) )

    # add peaks MS2
    pks_ms2 <- chromPeaksMS2( x = xdata, mtab_files = file_mtab[loop] )
    attr(xdata, "chromPeakMS2") <- pks_ms2
    xdata <- chromPeaksPostFilter(
      x = xdata, itab = itab, mtab = mtab, postfilter = postfilter, BPPARAM = BPPARAM 
      )

    saveRDS(xdata, file = file.path(dir_peak, f1b))
    file_XCMSnExp <- c( file_XCMSnExp, file.path(dir_peak, f1b) )
  }
  
  # flt <- list(postfilter = postfilter, noisefilter = noisefilter)
  # flt$x <- NULL
  # flt$postfilter <- NULL
  # flt$BPPARAM <- NULL
  
  list(itab = file_itab, 
       mtab = file_mtab, 
       XCMSnExp = file_XCMSnExp, 
       xcmsScanFilter = noisefilter, 
       peakPickingParam = param, 
       postFilter = postfilter)
}

#' define features for individual pipeline 
#' @param files the file paths of XCMSnExp
#' @param mtab_files the file path of meta tables, often returned by function \code{peakPicking}
#' @param pgParam the peak grouping parameter, passed to \code{groupChromPeaks}
#' @param rtParam the retention time adjustment parameters, passed to \code{adjustRtime}
#' @export
#' @importFrom stringr str_split_fixed
#' @importFrom S4Vectors DataFrame
#' @importFrom stats cutree dist hclust median
#' @return It returns an object of class \code{XCMSnExp}
#' 
defineFeatures <- function(files, mtab_files, rtParam = NULL, pgParam = PeakDensityParam(rep(1, length(files)))) {
  exps <- lapply(files, readRDS)
  xdata <- do.call(c, exps)
  xdata <- groupChromPeaks(xdata, param = pgParam)
  if (!is.null(rtParam) & length(files) > 1) {
    xdata <- adjustRtime(xdata, param = rtParam)
    xdata <- groupChromPeaks(xdata, param = pgParam)
  } 
  
  ##### redefine groups #####
  fd <- featureDefinitions(xdata)
  fd <- fd[, setdiff(colnames(fd), "ms_level"), drop = FALSE]
  fd <- as.data.frame(fd)
  cpd <- as.data.frame(chromPeaks(xdata))
  cpms2 <- lapply(exps, attr, "chromPeakMS2")
  cpms2 <- do.call(rbind, cpms2)
  if (nrow(cpms2) > 0)
    cpd <- cbind(cpd, cpms2)
  cpd$idx <- 1:nrow(cpd)
  chromPeaks(xdata) <- cpd
  
  writeFeature <- function(x, what) {
    df <- data.frame(
      mzmed = median(x$mz),
      mzmin = min(x$mz),
      mzmax = max(x$mz),
      rtmed = median(x$rt),
      rtmin = min(x$rt),
      rtmax = max(x$rt),
      npeaks = nrow(x),
      stringsAsFactors = FALSE
    )
    df$peakidx <- list(x$idx)
    if (any(!what %in% colnames(df)))
      stop(sprintf(
        "writeFeature: columns (%s) not available!", paste(setdiff(what, colnames(df)), collapse = ",")
      ))
    df[what]
  }
  
  splitPeaks <- function(x, h) {
    x <- hclust(dist(x, method = "man"), method = "single")
    cutree(x, h = h)
  }
  
  tt <- lapply(1:nrow(fd), function(i) {

    cc <- c("mzmed", "mzmin" ,"mzmax" ,"rtmed", "rtmin" , "rtmax", "npeaks", "peakidx")
    ii <- fd$peakidx[[i]]
    if (length(ii) == 1)
      return(fd[i, cc]) 
    
    tab <- cpd[ii, ]
    nrange <- 3
    splitMore <- TRUE
    while (splitMore) {
      sl <- splitPeaks(
        tab[, c("rt", "rtmin", "rtmax")], h = nrange*median(tab$rtmax - tab$rtmin)
      )
      tb <- table(sl, tab$sample) > 1
      if (sum(tb)/length(tb) > 0.25) {
        nrange <- nrange * 0.95
        splitMore <- TRUE
      } else
        splitMore  <- FALSE
    }    
    
    if (length(unique(sl)) == 1)
      return(fd[i, cc]) 
    
    st <- split(tab, sl)    
    p <- lapply(st, writeFeature, what = cc)
    do.call(rbind, p)
  })
  tq <- do.call(rbind, tt)
  n <- nrow(tq)
  rownames(tq) <- xcms:::.featureIDs(nrow(tq), prefix = "FT")
  tq <- DataFrame(tq, check.names = FALSE)
  class(tq$peakidx) <- mode(tq$peakidx)
  featureDefinitions(xdata) <- tq  
  ##### redefine groups end #####
  
  if (hasAdjustedRtime( xdata )) {
    rt <- rtime(xdata)
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
  xdata
}

#' Process raw file to generate the prunedXcmsSet
#' @param files mzXML/mzML files
#' @param peakPickingParam methods passed to \code{\link[xcms]{findChromPeaks}}
#' @param tmpdir directory for temporary files
#' @param postfilter post filter params to filter identified peaks, passed to \code{\link{chromPeaksPostFilter}}
#' @param noisefilter passed to \code{getMetaIntensityTable}
#' @param keepMS1 logical value. If the MS1 intensities should be kept.
#' @param keepEIC logical value. If the extract ion chromatogram (EIC) should be kept.
#' @param pheno phenotype data. It should be a data.frame with at least one column named as "file" to 
#'   list all the files passed to the function
#' @param RTAdjustParam methods passed to \code{\link[xcms]{adjustRtime}}
#' @param peakGroupParam methods passed to \code{\link[xcms]{groupChromPeaks}} for peak grouping
#' @param fillMissing logical; whether fill the missing values \code{\link[xcms]{fillChromPeaks}} 
#' @param mode Ionization mode, either 'pos' (default) or "neg"
#' @param ref MS1 and MS2 annotation data
#' @param ppmtol the mass tolerence given by parts per million (PPM)
#' @param mclapplyParam the parallel function, could be \code{mclapply} or \code{bplapply}, used in peak picking and annotation
#' @param bplapplyParam the parameters passed to biocparallel biolapply
#' @param peakPickingObj object save after the peak picking step, usually in "04_object_PeakPicking.RDS". If this is give the peak 
#'   picking step will be ignored.
#' @import parallel
#' @import BiocParallel
#' @import xcms
#' @export
#' @return An object of class "prunedXcmsSet"
#' 
runPrunedXcmsSet <- function( 
  
  # peakPicking
  files, 
  peakPickingParam = MatchedFilterParam(fwhm = 7.5), 
  tmpdir="xcmsViewerTemp",  
  postfilter = c(5, 1000), 
  noisefilter = c(
    ms1.noise = 100, ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 20,
    ms2.noise = 30, ms2.maxPeaks = 100, ms2.maxIdenticalInt = 6, 
    BPPARAM=bpparam()),

  keepMS1 = TRUE,
  keepEIC = TRUE,
  
  # phenotype data
  pheno = NULL,
  
  # feature identification
  RTAdjustParam = PeakGroupsParam(
    minFraction = 0.25,
    extraPeaks = 1,
    smooth = "loess",
    span = 0.5,
    family = "gaussian",
    peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
    subset = integer(),
    subsetAdjust = c("average", "previous")), 
  peakGroupParam = PeakDensityParam(seq_len(length(files))),
  
  fillMissing = FALSE,
  
  # feature annotation
  mode = c('pos', "neg")[1], 
  ref = NULL,
  ppmtol = 25,
  
  # parallel
  mclapplyParam = list(fun_parallel = mclapply, mc.cores = 1),
  bplapplyParam = MulticoreParam(workers = 1),
  peakPickingObj = NULL
) {
  
  if (!is.null(pheno))
    if (length(files) != nrow(pheno)) 
      stop("length(files) != nrow(pheno)!")
  
  # peakPicking and filter
  if (is.null(peakPickingObj)) {
    pp <- peakPicking( 
      files = files, param = peakPickingParam, tmpdir=tmpdir, 
      noisefilter = noisefilter , postfilter = postfilter, BPPARAM = bplapplyParam
      )
    
    saveRDS(pp, file = file.path(tmpdir, "04_object_PeakPicking.RDS"))
  } else
    pp <- peakPickingObj  
  
  # feature identification, the RT will also be identified
  df <- defineFeatures(
    files = pp$XCMSnExp, mtab_files = pp$mtab, rtParam = RTAdjustParam,  pgParam = peakGroupParam
  )
  
  # filling
  if (fillMissing) {
    df <- fillChromPeaks(df, msLevel = 1L, BPPARAM = bplapplyParam)
  } 

  cat("Extracting extended chrom peaks ...\n")
  peaks <- chromPeaks(df)
  peaks$ID <- rownames(peaks)
  peaks$is_filled <- chromPeakData(df)$is_filled
  peaks$ms_level <- chromPeakData(df)$ms_level
  cat(" done! \n")
  
  cat("Extracting extended feature table ...\n")
  fd <- as.data.frame(featureDefinitions(df))
  fd$ID <- rownames(fd)
  fd$annot_ms1 <- character(nrow(fd))
  fd$annot_ms2 <- character(nrow(fd))

  ########
  cat("Extracting Scans and EICs ...\n")
  itable <- c()
  mtable <- c()
  eicList <- list()
  for (i in 1:length(pp$itab)) {
    print(paste(i, length(pp$itab), sep = '/'))
    itab1 <- readRDS(pp$itab[i])
    mtab1 <- readRDS(pp$mtab[i])

    # === EIC start ===
    if (keepEIC) {
      scanList <- list( scanMeta = mtab1, scanIntensity = itab1 )      
      eicList[[i]] <- mclapplyParam$fun_parallel(1:nrow(fd), function(i) {
          f <- fd[i, ]
          mzoff <- 20 * 1e-06 * f$mzmed
          getEIC2(scanList, rt = c(f$rtmin - 30, f$rtmax + 30), 
                mz = c(f$mzmin - mzoff, f$mzmax + mzoff), na.rm = TRUE)
        }, mc.cores = mclapplyParam$mc.cores)
      names(eicList[[i]]) <- fd$ID
    } 
    # === EIC done ===

    if (!keepMS1) {    
      mtab1 <- mtab1[which(mtab1$msLevel > 1), ]
      itab1 <- itab1[itab1$ID %in% mtab1$ID, ]    
    }
    itable <- rbind(itable, itab1)
    mtable <- rbind(mtable, mtab1)
  }
  if (keepEIC) {
    eicList <- lapply(fd$ID, function(x) {
        do.call(rbind, lapply(eicList, "[[", x))
      })
    names(eicList) <- fd$ID
  }
  obj_xcmsScan <- new("xcmsScan", meta = mtable, intensity = itable, filter = pp$xcmsScanFilter, keepMS1 = TRUE)  
  
  ##### add consensus spectra #####
  cat("Extracting consensus spectra ...\n")
  cp <- mclapplyParam$fun_parallel( 1:nrow(fd), function(i, feat, peaks) {
    if (i %% 100 == 0)
      print(paste(i, nrow(fd)), sep = "/")
    pk <- peaks[feat$peakidx[[i]], ]
    pk <- pk[!pk$is_filled, ]
    unlist(strsplit(pk$ms2Scan, ";"))
    v <- obj_xcmsScan@intensity$ID %fin% unlist(strsplit(pk$ms2Scan, ";"))
    v <- obj_xcmsScan@intensity[v, ]
    if (nrow(v) == 0)
      return(c("", "0"))
    r <- consensusSpectrumLite(v, df = FALSE)
    c(r, attr(r, 'purity'))
  }, feat = fd, peaks = peaks, mc.cores = mclapplyParam$mc.cores )
  cp <- do.call(cbind, cp)
  
  fd <- cbind(fd, data.frame(
    ms2spectrum = cp[1, ],
    purity = as.numeric(cp[2, ]),
    stringsAsFactors = FALSE
  ))
  
  fd <- fd[, .xcmsViewerInternalObjects()$xcmsFeatureSet_fdata_name]
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
   
  cat("Summarizing experiment ...\n")
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
  if (is.numeric(peaks$validMS2))
    peaks$validMS2 <- as.logical(peaks$validMS2)
  obj_xcmsPeak <- new("xcmsPeak", 
                      table = peaks[, .xcmsViewerInternalObjects()$xcmsPeak_table_column],
                      param = pp$peakPickingParam,
                      postFilter = pp$postFilter)
  
  obj_xcmsFeatureSet <- new(
    "xcmsFeatureSet", 
    intensity = int, masterPeak = masterPeak, 
    phenoData = AnnotatedDataFrame(pheno),
    featureData = AnnotatedDataFrame(fd), 
    RTAdjustParam = RTAdjustParam, 
    peakGroupParam = peakGroupParam
  )

  res <- new("prunedXcmsSet", 
             featureSet = obj_xcmsFeatureSet,
             peak = obj_xcmsPeak, 
             scan = obj_xcmsScan, 
             EIC = eicList)
  saveRDS(res, file = file.path(tmpdir, "05_object_noAnnot.RDS"))

  ref <- prepRef(ref, mode = mode)
  amParam <- list(
    object = res, mode = mode, ref = ref, ppmtol = ppmtol
  )
  amParam <- c(amParam, mclapplyParam)
  saveRDS(amParam, file = file.path(tmpdir, "06_object_annot_param.RDS"))
  do.call(annotateMetabolite, amParam)
}

