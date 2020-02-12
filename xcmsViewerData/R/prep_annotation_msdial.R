#' Processing MSP files of msDial to list of meta table for ms1 and ms2 and peak list
#' @param mspFile .msp file
#' @param db the database prefix, used for creating IDs
#' @param predicted predicted spectra or expermental
#' @param cores passed to mclapply
#' @importFrom parallel mclapply
#' @importFrom stringr str_split_fixed
#' @import fastmatch
#' @export
#' 
#' 
prep_annotation_msdial <- function(mspFile, db, predicted = c("false", "true")[1], cores = 1) {
  
  msp <- readLines(mspFile)
  msp <- split(msp, cumsum(msp==""))
  
  mspl <- mclapply(msp, mc.cores = cores, function(msp1) {
    
    l <- list(meta_ms1 = NULL, meta_ms1 = NULL, peak = NULL)
    if (length(pr <- grep("^Num\ Peaks", msp1)) == 0)
      return(l)
    ir <- 1:pr
    meta <- stringr::str_split_fixed(setdiff(msp1[ir], ""), pattern = ": ", 2)
    meta <- structure(trimws(meta[, 2]), names = trimws(meta[, 1]))
    
    inchi <- paste0("INCHIKEY:", .getVal(meta, "INCHIKEY"))
    
    ## ms1 meta
    l$meta_ms1 <- c("name" = .getVal(meta, "NAME"), 
                    "chemical_formula" = .getVal(meta, "FORMULA"), 
                    "average_molecular_weight" = NA, 
                    "monoisotopic_molecular_weight" = NA, 
                    "cas_registry_number" = inchi, 
                    "accession" = NA, 
                    prec_type = .getVal(meta, "PRECURSORTYPE"),
                    prec_mass = .getVal(meta, "PRECURSORMZ")
                    )
    
    ## msms meta
    l$meta_ms2 <- c("id"=NA, "notes" = .getVal(meta, "Comment"), 
                    "sample-concentration" = NA, "solvent" = NA, 
                    "sample-mass" = NA,
                    "sample-assessment" = NA, "spectra-assessment" = NA, 
                    "sample-source" = NA, "collection-date" = NA, 
                    "instrument-type" = .getVal(meta, "INSTRUMENTTYPE"),
                    "peak-counter" = .getVal(meta, "Num Peaks"), "created-at" = NA, "updated-at" = NA,
                    "mono-mass" = NA, "energy-field" = .getVal(meta, "COLLISIONENERGY"),
                    "collision-energy-level" = .getVal(meta, "COLLISIONENERGY"), 
                    "collision-energy-voltage" = .getVal(meta, "COLLISIONENERGY"),
                    "ionization-mode" = .getVal(meta, "IONMODE"),
                    "sample-concentration-units" = NA, "sample-mass-units" = NA,
                    "predicted" = predicted, "structure-id" = inchi, "splash-key" = NA)
    
    ## msms peak
    if (!is.na(.getVal(meta, "MSLEVEL")))
      if (.getVal(meta, "MSLEVEL") == "MS1")
        l$peak <- NULL
    int <- stringr::str_split_fixed(msp1[-ir], pattern = "\t", 3)
    int <- data.frame("mz" = int[, 1], "intensity" = int[, 2], stringsAsFactors = FALSE)
    l$peak <- stats::na.omit(int)
    l
  })
  mspl <- mspl[!sapply(mspl, function(x) is.null(x$meta_ms1))]
  
  meta_ms1 <- data.frame(t(sapply(mspl, "[[", "meta_ms1")), stringsAsFactors = FALSE)
  meta_ms1$average_molecular_weight <- as.numeric(meta_ms1$average_molecular_weight)
  meta_ms1$monoisotopic_molecular_weight <- as.numeric(meta_ms1$monoisotopic_molecular_weight)
  meta_ms1$prec_mass <- as.numeric(meta_ms1$prec_mass)
  meta_ms1$monoisotopic_molecular_weight <- .toMonoisoMass(meta_ms1$prec_type, meta_ms1$prec_mass)
  meta_ms1$prec_type <- NULL
  meta_ms1$prec_mass <- NULL
  
  meta_ms2 <- data.frame(t(sapply(mspl, "[[", "meta_ms2")), stringsAsFactors = FALSE, check.names = FALSE)
  meta_ms2$id <- as.character(1:nrow(meta_ms2))
  meta_ms2[["database-id"]] <- .msms_id_creater(num = meta_ms2$id, db = db)
  meta_ms2$id2 <- meta_ms2[["database-id"]]
  meta_ms2$"sample-concentration" <- suppressWarnings(as.numeric(meta_ms2$"sample-concentration"))
  meta_ms2$"sample-mass" <- suppressWarnings(as.numeric(meta_ms2$"sample-mass"))
  meta_ms2$"mono-mass" <- suppressWarnings(as.numeric(meta_ms2$"mono-mass"))
  meta_ms2$"collision-energy-voltage" <- suppressWarnings(as.numeric(meta_ms2$"collision-energy-voltage"))
  
  meta_ms1$accession <- meta_ms2$id2
  meta_ms1 <- meta_ms1[!is.na(meta_ms1$monoisotopic_molecular_weight), ]
  
  peaks <- lapply(mspl, function(x)  {
    x <- x$peak
    if (length(x) == 0)
      return(NULL)
    if (nrow(x) == 0)
      return(NULL)
    data.frame(
      id = factor(NA), intensity = as.numeric(x$intensity),
      mz = as.numeric(x$mz), molecule.id = factor(NA), ms.ms.id = factor(NA), 
      stringsAsFactors = FALSE)
  })
  names(peaks) <- meta_ms2$id2
  
  ii <- !sapply(peaks, is.null) && names(peaks) %fin% meta_ms1$accession 
  list(meta_ms1 = meta_ms1,
       meta_ms2 = meta_ms2[ii, ],
       peaks = peaks[ii]
  )
}

