#' Prepare hmdb data for mass annotation in xcms file
#' @param x an xml file downloaded from hmdb database or an RDS file reading the RDS.
#' @importFrom xml2 read_xml as_list
#' @export
#' 

prep_annotation_hmdb <- function(x) {
  
  if (grepl(".xml$", x, ignore.case = TRUE)) {
    cat("Reading xml file ...")
    xf <- read_xml(x)
    xflist <- as_list(xf)
  } else if (grepl(".rds$", x, ignore.case = TRUE)) {
    cat("Reading RDS file ...")
    xflist <- readRDS(x)
  } else {
    stop("Unknown file format for x")
  }
  
  ext <- .ms_meta_header()
  ss <- sapply(xflist$hmdb, function(met) {
    sapply(ext, function(i) {
      if (!(i %in% names(met)))
        return(NA)
      paste(unlist(met[[i]]), collapse = ";")
    })
  })
  ss <- t(ss)
  ss <- apply(ss, 2, trimws)
  ext[ext == "monisotopic_molecular_weight"] <- "monoisotopic_molecular_weight"
  colnames(ss) <- ext  
  ss <- data.frame(ss, stringsAsFactors = FALSE)
  
  ss$average_molecular_weight <- as.numeric(ss$average_molecular_weight)
  ss$monoisotopic_molecular_weight <- as.numeric(ss$average_molecular_weight) 
  ss
}