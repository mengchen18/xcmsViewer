#' parse rt and mz range used by shiny app
#' @param x a character to be parsed, such as "123-133" or "133"
#' @param expand when x is parsed to a single value, how to expand the range
parseRange <- function(x, expand = 0.05) {
  v <- as.numeric(strsplit(x, "-")[[1]])
  if (any(is.na(v))) {
    warning('Unknown format x!')
    v <- c(0, 0)
  } else if (length(v) == 0) {
    v <- c(-Inf, Inf)
  } else if (length(v) == 2) {
    if (v[1] >= v[2]) {
      warning("Wrong range definition!")
      v <- c(0, 0)
    }
  } else if (length(v) == 1) {
    v <- c(max(0, v-expand), v+expand)
  } else {
    warning('Unknown format x!')
    v <- c(0, 0)
  }
  v
}

#' format feature or peak table
#' @param tab table
formatFeaturePeakTab <- function(tab) {
  numc <- which(sapply(tab, function(x) is.numeric(x) && !is.integer(x)))
  dt <- DT::datatable( 
    tab,
    selection =  list(mode="single", selected = 1),
    rownames = FALSE,
    filter = "top",
    class="table-bordered compact",
    caption = "Associated peaks",
    options = list(scrollX = TRUE, scrollY = "180px", paging = FALSE, dom = 't')
  )
  dt <- DT::formatRound(dt, columns = numc, digits = 3)
  DT::formatStyle(dt, columns = 1:ncol(tab), fontSize = '85%')
}

#' prepare MGF file to be download
#' @param scan scan list, a data.frame with three columns - mz, intensity, ID
#' @param scanMeta scan meta data
#' @param cons consensusMS2Spectrum - a data.frame with two columns, mz and intensity
#' @param mode 1+ or 1-
#' @param featureId ID
#' @examples 
#' # x is obj() from xcmsAnnotationTab_module
#' # mgf <- prepareMGF(
#' # scan = x$ms2scan,
#' # scanMeta = x$ms2scanMeta,
#' # cons = x$consensusMS2Spectrum,
#' # mode = "pos",
#' # featureId = x$featureid
#' # )
#' # writeLines(mgf, con = "spectrum.mfg")
#' 
prepareMGF <- function(scan, scanMeta, cons, mode, featureId) {
  
  scanMeta <- scanMeta[scanMeta$ID %in% scan$ID, ]
  if (nrow(scanMeta) == 0)
    return(NULL)
  
  peaks <- paste(cons$mz, cons$intensity)
  meta <- sprintf(
    "Feature:%s|RT:%s|TIC:%s|ID:Consensus_MS2_spectrum_of_%s", 
    featureId,
    median(scanMeta$rt, na.rm = TRUE), 
    sum(scanMeta$tic, na.rm = TRUE),
    nrow(scanMeta))
  
  pc <- c(
    "BEGIN IONS",
    paste0("PEPMASS=", median(scanMeta$precMz, na.rm = TRUE)),
    paste0("CHARGE=", mode),
    "MSLEVEL=2",
    paste0("TITLE=", meta),
    peaks,
    "END IONS",
    "")
  
  pks <- lapply(1:nrow(scanMeta), function(i) {
    m1 <- scanMeta[i, ]
    peaks <- scan[scan$ID == m1$ID, ]
    if (nrow(peaks) == 0)
      return(NULL)
    peaks <- paste(peaks$mz, peaks$intensity)
    
    meta <- sprintf("Feature:%s|RT:%s|TIC:%s|ID:%s", featureId, m1$rt, m1$tic, m1$ID)
    c("BEGIN IONS",
      paste0("PEPMASS=", m1$precMz),
      paste0("CHARGE=", mode),
      "MSLEVEL=2",
      paste0("TITLE=", meta),
      peaks,
      "END IONS",
      "")
  })
  
  c(pc, unlist(pks))
}


#' Convert spectra from string format to data.frame format
#' @param x a character (vector) of spectra in string format,
#'   such as "78.454;123.333;332.232|50;1500;400"
#' @return a data.frame with two columns "mz" and "intensity"
#' 
str2spectra <- function(x) {
  
  if (any(is.na(x))) {
    warning("NA input in str2spectra!")
    return(NULL)
  }
  .l <- strsplit(x, "\\|")
  
  .l <- lapply(.l, function(xx) {
    
    if (any(is.na(xx))) {
      warning("NA input in str2spectra!")
      return(NULL)
    }
    
    if (length(xx) == 0)
      return(NULL)
    
    ls <- strsplit(xx, ";")
    
    if (length(ls[[1]]) != length(ls[[2]])) {
      warning("Length of spectrum intensity is not the same as spectrum mass!")
      return(NULL)
    }
    
    r <- NULL
    if (length(ls) > 0) {
      r <- data.frame(
        mz = as.numeric(ls[[1]]),
        intensity = as.numeric(ls[[2]])
      )
    }
    
    r
  })
  
  if (length(.l) == 1)
    .l <- .l[[1]]
  .l
}

