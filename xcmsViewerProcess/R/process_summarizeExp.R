#' Summarize XCMSnExp data
#' @param x an object of XCMSnExp
#' @param mode mode. pos or neg
#' @param massTab a data.frame of mass table having at least one column named as "monoisotopic_molecular_weight"
#' @param refSpectra reference MS2 spectra used annotate the experimental MS2 spectra
#' @param ... other parameters passed to bplapply 
#' @importFrom fastmatch %fin%
#' @import MAIT
#' @import xcms
#' @import BiocParallel
#' @export
#' 
summarizeExp <- function(x, mode = c("pos", "neg")[1], massTab=NULL, refSpectra = NULL, ...) {
  
  mtable <- asMetaTable(x)

  cat("Extracting peak intensity table and meta table ...")
  itable <- asIntensityTable(x, ...)  
  cat(" done! \n")
  
  cat("Extracting extended chrom peaks ...\n")
  peaks <- extendedChromPeaks(x, mtab=mtable, itab=itable, ...)
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
  peaks <- peaks[, c("ID",  "QC", "sample", "masterPeak", "into", 
    "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "intb", "maxo", "sn", 
    "ms2Scan", "rsq", "rtgap", "intgap", "rtintgap", "truncated", "b"), ]
  cat(" done! \n")
  
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
  
  if (!is.null(massTab)) {
    ll <- bplapply(1:nrow(features$meta), function(i) {
      massQuery(m = features$meta$mzmed[i], tolppm = 10, refTab = massTab, 
        addTable = at, ID = rownames(features$meta)[i]) 
    }, ...)
    names(ll) <- features$meta$ID
  } else 
    ll <- NULL
  gc()
  cat("done!\n")
  
  ail <- NULL 
  ms2a  <- NULL
  if (!is.null( refSpectra )) {
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
      ail <- NULL 
      ms2a  <- NULL
    } else {
      ms2a <- do.call(rbind, ms2an)
      colnames(ms2a) <- gsub("-", "_", colnames(ms2a))
      ap <- unique(ms2a$annot_peaks)

      fmatch(ap, refSpectra$meta$id2)

      ail <- list(
        # meta = refSpectra$meta[fmatch(ap, refSpectra$meta$id2), ],
        meta = refSpectra$meta[refSpectra$meta$id2 %fin% ap, ],
        peakList = refSpectra$peakList[ap]
      )
      colnames(ail$meta) <- gsub("-", "_", colnames(ail$meta))
    }
    gc()
    cat("done!\n")
    
    cat("Combining MS2 and MS1 annotations ... \n")
    vx <- lapply(features$meta$ID, function(id) {
      ## annotation by mass comparison
      s_massAt <- ll[[id]]
      i <- fmatch(id, features$meta$ID)
      if (is.null(s_massAt))
        return(NULL)
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
    cat("Done!\n")
  }
  
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