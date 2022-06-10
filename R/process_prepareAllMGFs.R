#' prepare MGF file to be download
#' @param obj an object of prunedXcmsSet, returned by prepVeiwerData
#' @param file output file
#' @export

prepareAllMGFs <- function(obj, file) {
  cat("Extracting consensus spectra ...\n")
  fd <- getFeatureMeta(obj)
  cspec <- str2spectra(fd$`General|Extended|ms2spectrum`)
  names(cspec) <- fd$`General|All|ID`

  cat("Extracting scan info ...\n")
  s_meta <- getScanMeta(obj) 
  s_int <- getScan2(obj) 
  s_meta <- s_meta[s_meta$ID %in% s_int$ID, ]
  s_int <- split(s_int, s_int$ID)
  s_int <- s_int[sapply(s_int, nrow) > 0]
  
  cat("Splitting ...\n")
  ids <- sapply(fd$`General|All|ID`, function(x) getScanIDFromFeatureID(object = obj, x))
  sin <- names(s_int)
  ids <- sapply(ids, function(x) x[x %fin% sin])
  s_meta_list <- lapply(ids, function(x) s_meta[fmatch(x, s_meta$ID), ])
  
  mod <- getIonMode(obj)
  if (grepl("pos", mod, ignore.case = TRUE)) (
    md1 <- "1+"
  ) else if (grepl("neg", mod, ignore.case = TRUE)) {
    md1 <- "1-"
  } else
    stop("Unknown mode, prepare MFGs!")    
  
  cspec <- cspec[!sapply(cspec, is.null)]
  if (length(cspec) == 0)
    return(NULL)

  cat("Creating MGFs ...\n")
  ll <- lapply(names(cspec), function(ft) {
    m1 <- s_meta_list[[ft]]
    m1s <- do.call(rbind, s_int[m1$ID])
    prepareMGF(scan = m1s, scanMeta = m1, cons = cspec[[ft]], featureId = ft, mode = md1)
  })
  ll <- do.call(c, ll)
  
  cat("Writing ...\n")
  writeLines(ll, con = file)
  
  invisible(ll)
}