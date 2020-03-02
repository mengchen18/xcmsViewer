#' @describeIn summarizeExp summarize XCMSnExp to simple table/list format
#' @export

xcms_summarize <- function(
  x, 
  ms1.noise = 100, ms1.maxPeaks = Inf, ms1.maxIdenticalInt = 20,
  ms2.noise = 30, ms2.maxPeaks = 100, ms2.maxIdenticalInt = 6, 
  QC = TRUE, itable = NULL, mtable = NULL, ...) {

  if (is.null(mtable))
    mtable <- asMetaTable(x)
  
  if (is.null(itable)) {
  cat("Extracting peak intensity table table ...")
  itable <- asIntensityTable(
    x,  
    ms1.noise = ms1.noise, ms1.maxPeaks = ms1.maxPeaks, ms1.maxIdenticalInt = ms1.maxIdenticalInt,
    ms2.noise = ms2.noise, ms2.maxPeaks = ms2.maxPeaks, ms2.maxIdenticalInt = ms2.maxIdenticalInt, ...)  
  cat(" done! \n")
  } else 
    message("itable is given, ms1 and ms2 parameters for itable are not used!")
  
  cat("Extracting extended chrom peaks ...\n")
  gc()
  peaks <- extendedChromPeaks(x, mtab=mtable, itab=itable, QC = QC, ...)
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
    scanMetaTab = mtable,
    scanIntensityTab = itable,
    pheno = pheno
  )
  res
}