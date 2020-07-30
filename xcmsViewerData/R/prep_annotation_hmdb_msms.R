#' parse MSMMS XML files downloaded from HMDB
#' @param xmlFolder the xml folder
#' @param ... other parameters passed to bplapply
#' @importFrom flatxml fxml_importXMLFlat
#' @importFrom stringr str_pad
#' @import BiocParallel
#' @export

prep_annotation_hmdb_msms <- function(xmlFolder, ...) {
  
  metaHeader <- .msms_meta_header()
  
  xmls <- list.files(xmlFolder, pattern = ".xml", full.names = TRUE)
  lx <- length(xmls)
  fl <- bplapply(1:length(xmls), function(i) {
    if (i%%10000 == 0)
      cat(sprintf("%s out of %s has finished!\n", i, lx))
    f1 <- xmls[i]
    ff <- fxml_importXMLFlat(f1)
    ff <- ff[!is.na(ff$value.), ]
    
    i <- is.na(ff$level3) & is.na(ff$level4)
    meta <- structure(ff$value.[i], names = ff$elem.[i])
    meta[metaHeader]
    
    i <- ff$level3 == "ms-ms-peak" & ff$level2 == "ms-ms-peaks"
    peakList <- split(ff$value.[i], ff$level4[i])
    if (length(unique(sapply(peakList, length))) > 1)
      stop("peak list has different length.")
    peakList <- as.data.frame(peakList, stringsAsFactor = FALSE)
    peakList$intensity <- as.numeric(as.character(peakList$intensity))
    peakList$mass.charge <- as.numeric(as.character(peakList$mass.charge))
    colnames(peakList)[colnames(peakList) == "mass.charge"] <- "mz"
    list(meta = meta, peakList = peakList)
  }, ...)
  
  meta <- t(sapply(fl, "[[", "meta"))
  meta[meta == "true"] <- NA
  meta <- as.data.frame(meta, stringsAsFactor = FALSE)
  for (nn in names(meta))
    meta[[nn]] <- as.character(meta[[nn]])
  numcol <- c("sample-concentration", "sample-mass", "mono-mass", "collision-energy-voltage")
  meta[numcol] <- lapply(meta[numcol], function(x) as.numeric(as.character(x)))

  meta <- data.frame(
    "mono_mass" = meta$"mono-mass",
    "solvent" = meta$solvent,
    "collision_energy_level" = meta$"collision-energy-level",
    "collision_energy_voltage" = meta$"collision-energy-voltage",
    "ionization_mode" = meta$"ionization-mode",
    "predicted" = meta$"predicted",
    "database" = "HMDB",
    "database_id" = meta$"database-id",
    "structure_id" = meta$"structure-id",
    stringsAsFactor = FALSE
    )
  
  peakList <- lapply(fl, "[[", "peakList")
  names(peakList) <- meta$id2
  
  ir <- sapply(peakList, nrow) > 0
  meta <- meta[ir, ]
  peakList <- peakList[ir]
  
  list(meta = meta, peakList = peakList)
}


