#' Validation functions
#' @noRd
#' @importFrom utils data read.delim
.validate_df_general <- function(x, name, class, what = c("class", "empty", "col_class"), 
                                other_class=NULL, verbose = TRUE) {
  
  if (verbose)
    cat(sprintf("validating %s table ...", name))
  
  # check class and empty
  if ("class" %fin% what) {
    if (!inherits(x, "data.frame"))
      message("Possible problem: Table should be a data.frame!")
  }
  
  if ("empty" %fin% what) {
    if (nrow(x) == 0)
      message("Possible problem: Seems the peaks table is empty!")
  }
  
  if ("col_class" %fin% what) {
    cls <- class
    
    # check column names
    dffcol <- setdiff(cls$header, colnames(x))
    if (length(dffcol) > 0) {
      dffcol <- paste(dffcol, collapse = ", ")
      message(sprintf("Possible problem: These columns are missed in the %s table: %s", name, dffcol))
    }
    
    # check column class
    xcls <- sapply(x, class)[cls$header]
    i <- which(xcls != cls$class)
    if (length(i) > 0) {
      for (ii in i) {
        abn <- xcls[ii]
        message(
          sprintf("Error: Column %s should be an %s class, but not %s!", names(abn), cls$class[ii], abn)
        )
      }
      message("Possible problem: Check column class!")
    }
    
    dffcol_2 <- setdiff(colnames(x), cls$header) 
    
    if (length(dffcol_2) > 0) {
      if (is.null(other_class))
        message("Possible problem: Extra columns present but the their class is not defined in other class!")
      v <- unique(sapply(x[, dffcol_2], function(x) inherits(x, other_class)))
      if (any(!v))
        message(sprintf("Possible problem: Extra columns are not the class of %s.", other_class))
    }
  }
  
  if (verbose)
    cat("done!\n")
}

.validate_peaks <- function(x) {
  
  cls <- data.frame(
    header =  c(
      "ID", "QC", "sample", "masterPeak", "into", 
      "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", 
      "intb", "maxo", "sn", "ms2Scan", "rsq", 
      "rtgap", "intgap", "rtintgap", "truncated", "b"),
    class = c(
      "character", "character", "factor", "character", 
      "numeric", "numeric", "numeric", "numeric", 
      "numeric", "numeric", "numeric", "numeric", 
      "numeric", "numeric", "character", "numeric", 
      "numeric", "numeric", "numeric", "numeric", 
      "numeric"),
    stringsAsFactors = FALSE
  )
  
  .validate_df_general(x, "peaks", class = cls)
}


.validate_scanMetaTab <- function(x) {
  
  cls <- data.frame(
    header = c("scanNum", "acquisitionNum", "rt", "tic", "peakCount", "msLevel", 
               "fromFile", "precScanNum", "precMz", "precCharge", "precIntensity", "ID"),
    class = c("integer", "integer", "numeric", "numeric", "integer", "integer", 
              "integer", "integer", "numeric", "integer", "numeric", "character"),
    stringsAsFactors = FALSE
  )
  
  .validate_df_general(x, "scanMetaTab", class = cls)
  
}


.validate_scanIntensityTab <- function(x) {
  
  cls <- data.frame(
    header = c("mz", "intensity", "ID"),
    class = c("numeric", "numeric", "character"),
    stringsAsFactors = FALSE
  )
  
  .validate_df_general(x, "scanIntensityTab", class = cls)
}


.validate_annotationFragment <- function(x) {
  
  cls <- data.frame(
    header = c("annot_peaks", "query_peaks", "cos", "database_id"),
    class = c("character", "character", "numeric", "character" ),
    stringsAsFactors = FALSE
  )
  
  .validate_df_general(x, "annotationFragment", class = cls)
}

.validate_pheno <- function(x) {
  .validate_df_general(x, what = c("class", "empty"), name = "pheno")
}

.validate_features_meta <- function(x) {
  cls <- data.frame(
    header = c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", 
               "npeaks", "peakidx", "ID", "QC", "Annotation"), 
    class = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
              "numeric", "AsIs", "character", "character", "character"),
    stringsAsFactors = FALSE
  )
  .validate_df_general(x, name = "features$meta", class=cls, other_class="numeric")
}

.validate_features_intensities <- function(x) {
  cls <- data.frame(
    header = "ID", 
    class = "character",
    stringsAsFactors = FALSE
  )
  .validate_df_general(x, name = "features$intensities", class=cls, other_class="numeric")
}

.validate_features_masterPeaks <- function(x) {
  cls <- data.frame(
    header = "ID", 
    class = "character",
    stringsAsFactors = FALSE
  )
  .validate_df_general(x, name = "features$masterPeaks", class=cls, other_class="character")
}

.validate_annotationMass <- function(x) {
  cat('Validating annotationMass table ...')
  cls <- read.delim(textConnection(
    "
header class
name character
chemical_formula character
average_molecular_weight numeric
monoisotopic_molecular_weight numeric
cas_registry_number character
accession character
adduct character
massdiff numeric
massWithAdduct numeric
massQueried numeric
deltaPPM numeric
score numeric
ID character
annot_peaks character
query_peaks character
cos numeric
database_id character
atScore numeric
"), sep = " ", stringsAsFactors = FALSE)
  
  for (i in names(x)) {
    xx <- x[[i]]
    if (!is.null(xx)) {
      if (any(xx$ID != i))
        message("Possible problem: Validating annotationMassID doesn't map!")
      .validate_df_general(
        x=xx, name = paste0("annotationMass$", i), class=cls, verbose = FALSE)
    }
  }
  cat("Done!\n")
}

.validate_matchedRefFragments_meta <- function(x) {
  cls <- read.delim(textConnection(
    "
header class
id character
notes character
sample_concentration numeric
solvent character
sample_mass numeric
sample_assessment character
spectra_assessment character
sample_source character
collection_date character
instrument_type character
peak_counter character
created_at character
updated_at character
mono_mass numeric
energy_field character
collision_energy_level character
collision_energy_voltage numeric
ionization_mode character
sample_concentration_units character
sample_mass_units character
predicted character
structure_id character
splash_key character
database_id character
id2 character"
  ), sep = " ", stringsAsFactors = FALSE)
  
  .validate_df_general(x, name = "matchedRefFragments$meta", class=cls)
  
}



.validate_matchedRefFragments_peakList <- function(x) {
  
  cat("validating matchedRefFragments_peakList ...")
  cls <- data.frame(
    header = c("intensity", "mz"),
    class = c("numeric","numeric"),
    stringsAsFactors = FALSE
  )
  
  for (i in names(x)) {
    xx <- x[[i]]
    .validate_df_general(xx, name = paste0("matchedRefFragments$peakList$", i), 
                        class=cls, other_class="factor", verbose = FALSE)  
  }
  cat("done!\n")
  
}


.validate_features <- function(x) {
  
  it <- intersect(x$meta$ID, x$intensities$ID)
  it <- intersect(x$masterPeaks$ID, it)
  if (length(it) != nrow(x$meta))
    message("Possible problem: Validating features - ID columns in feature tables cannot be mapped!")
  
}


.validate_matchedRefFragments <- function(x) {
  if (any(!names(x$peakList) %fin% x$meta$id2))
    message("Possible problem: Validating matchedRefFragments - There are names in the peakList not recorded in the meta id2!")
}


.validate_all <- function(x) {
  
  if (!all(names(x$annotationMass) %fin% x$features$meta$ID) ||
      !all(x$features$meta$ID %fin% names(x$annotationMass)))
    message("Possible problem: Different elements: names of annotationMass and features$meta$ID")
  
  mstpk <- na.omit(unlist(x$features$masterPeaks[, setdiff(colnames(x$features$masterPeaks), "ID")]))
  if (any(!mstpk %fin% x$peaks$ID))
    message("Possible problem: Master peaks is not listed in the peaks table!")
  
  ###
  if (max(unlist(x$features$meta$peakidx)) != length(x$peaks$ID))
    message("Possible problem: features$meta$peakidx doens't match x$peaks$ID!")
  ###
  
  all_database_id <- unique(unlist(lapply(x$annotationMass, function(x) x$accession)))
  if (any(!x$annotationFragment$database_id %fin% all_database_id))
    message("Possible problem: Metabolite ID of MS2 fragments is not in the MS1 annotation (annotaitonMass) table!")
  
  if (any(!unique(unlist(strsplit(x$peaks$ms2Scan, ";"))) %fin% x$scanMetaTab$ID))
    message("Possible problem: peaks$ms2scan contains ID doesn't present in scanMetaTable$ID!")
  
  if (any(!x$scanIntensityTab$ID %fin% x$scanMetaTab$ID))
    message("Possible problem: scanIntensityTab$ID contains ID doesn't present in scanMetaTable$ID!")
  
  if (any(!x$annotationFragment$query_peaks %fin% x$scanMetaTab$ID))
    message("Possible problem: annotationFragment$query_peaks contains ID doesn't present in scanMetaTable$ID!")
  
  if (!all(x$annotationFragment$annot_peaks %fin% x$matchedRefFragments$meta$id2))
    message("Possible problem: annotationFragment$annot_peaks contains ID doesn't present in matchedRefFragments$meta$id2!")
  
  if (!all(x$annotationFragment$database_id %fin% x$matchedRefFragments$meta$database_id))
    message("Possible problem: annotationFragment$database_id contains ID doesn't present in matchedRefFragments$meta$database_id!")
}


.validate_statsXCMS <- function(x) {
  
  .validate_features_meta(x$features$meta)
  .validate_features_intensities(x$features$intensities)
  .validate_features_masterPeaks(x$features$masterPeaks)
  .validate_features(x$features)
  
  .validate_peaks(x$peaks)
  
  .validate_scanMetaTab(x$scanMetaTab)
  
  .validate_scanIntensityTab(x$scanIntensityTab)
  
  .validate_annotationMass(x$annotationMass) ####
  
  .validate_annotationFragment(x$annotationFragment)
  
  .validate_matchedRefFragments_meta(x$matchedRefFragments$meta)
  .validate_matchedRefFragments_peakList(x$matchedRefFragments$peakList)
  
  .validate_matchedRefFragments(x$matchedRefFragments)
  
  .validate_pheno(x$pheno)
  .validate_all(x)
  
}




