#' @import fastmatch, Biobase
#' 
############################# global variables ########################

.xcmsViewerInternalObjects <- new.env()

# xcmsFeatureSet
.xcmsViewerInternalObjects$xcmsFeatureSet_fdata_name <- c(
  "mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", "npeaks", 
  "peakidx", "ID", "annot_ms1", "annot_ms2")
.xcmsViewerInternalObjects$xcmsFeatureSet_fdata_name_class <- c(
  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
  "numeric", "AsIs", "character", "character", "character")

# xcmsScan
.xcmsViewerInternalObjects$xcmsScan_meta_column <- c(
  "scanNum", "acquisitionNum", "rt", "tic", "peakCount", "msLevel", "fromFile", 
  "precScanNum", "precMz", "precCharge", "precIntensity", "ID")
.xcmsViewerInternalObjects$xcmsScan_meta_column_class <- c(
  "integer", "integer", "numeric", "numeric", "integer", "integer", "integer", 
  "integer", "numeric", "integer", "numeric", "character")
.xcmsViewerInternalObjects$xcmsScan_intensity_column <- c("mz", "intensity", "ID")
.xcmsViewerInternalObjects$xcmsScan_intensity_column_class <- c("numeric", "numeric", "character")

# xcmsPeak
.xcmsViewerInternalObjects$xcmsPeak_table_column <- c(
  "ID", "sample", "masterPeak", "into", "mz", "mzmin", 
  "mzmax", "rt", "rtmin", "rtmax", "sn", "ms2Scan")
.xcmsViewerInternalObjects$xcmsPeak_table_column_class <- c(
  "character", "integer", "character", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", "numeric", "character")


############################ prunedXcmsSet ############################
setClass(
  "prunedXcmsSet", 
  slot = c(
    featureSet = "xcmsFeatureSet",
    peak = "xcmsPeak",
    scan = "xcmsScan"
  )
)


############################ xcmsFeatureSet ############################

setClass(
  "xcmsFeatureSet", 
  contains = "ExpressionSet"
)

setMethod("initialize", "xcmsFeatureSet", function(.Object, intensity, masterPeak, ...) {
  
  mv <- new.env()
  if (missing(intensity))
    mv$exprs <- matrix(numeric(0), 0, 0) else
      mv$exprs <- intensity
    if (missing(masterPeak))
      mv$masterPeak <- matrix(character(0), 0, 0) else
        mv$masterPeak <- masterPeak
      
      obj <- new("ExpressionSet", assayData = mv, ...)
      ll <- list(...)
      
      if (!"featureData" %in% names(ll)) {
        nr <- nrow(exprs(obj))
        fData(obj) <- data.frame(
          mzmed = numeric(nr),
          mzmin = numeric(nr),
          mzmax = numeric(nr),
          rtmed = numeric(nr),
          rtmin = numeric(nr),
          rtmax = numeric(nr),
          npeaks = numeric(nr),
          peakidx = I(rep(NA, nr)),
          ID = character(nr),
          annot_ms1 = character(nr),
          annot_ms2 = character(nr),
          stringsAsFactors = FALSE)
      }
      
      class(obj) <- "xcmsFeatureSet"
      validObject(obj)
      obj
})

setValidity("xcmsFeatureSet", function(object) {
  pd <- fData(object)
  .validRestrictedDataFrame(
    pd,
    name = .xcmsViewerInternalObjects$xcmsFeatureSet_fdata_name,
    class = .xcmsViewerInternalObjects$xcmsFeatureSet_fdata_name_class,
    str = "xcmsFeatureSet@fData",
    contain = TRUE                            )
})



############################ xcmsPeak ############################

setClassUnion(
  "chromPeakParamORnull", 
  c("NULL",
    "CentWaveParam",
    "CentWavePredIsoParam",
    "MatchedFilterParam",
    "MassifquantParam",
    "MSWParam"))
setClassUnion("numericORnull", c("NULL", "numeric"))

setClass(
  "xcmsPeak", 
  slots = c(
    table = "data.frame",
    param = "chromPeakParamORnull",
    postFilter = "numericORnull"
  )
)

setMethod("initialize", "xcmsPeak", function(.Object, table, param, postFilter) {
  
  if (missing(table)) {
    .Object@table = data.frame(
      "ID" = character(0),
      "sample" = integer(0),
      "masterPeak" = character(0),
      "into" = numeric(0),
      "mz" = numeric(0),
      "mzmin" = numeric(0),
      "mzmax" = numeric(0),
      "rt" = numeric(0),
      "rtmin" = numeric(0),
      "rtmax" = numeric(0),
      "sn" = numeric(0),
      "ms2Scan" = character(0),
      stringsAsFactors = FALSE
    )
  } else 
    .Object@table <- table[.xcmsViewerInternalObjects$xcmsPeak_table_column]
  
  if (missing(param))
    .Object@param <- NULL else
      .Object@param <- param
    
    if (missing(postFilter))
      .Object@postFilter <- NULL else
        .Object@postFilter <- postFilter
      
      validObject(.Object)
      return(.Object)
      
})

setValidity("xcmsPeak", function(object) {
  .validRestrictedDataFrame(
    df = object@table, 
    name = .xcmsViewerInternalObjects$xcmsPeak_table_column, 
    class = .xcmsViewerInternalObjects$xcmsPeak_table_column_class,
    str = "object@table"
  )
})


############################ xcmsScan ############################

setClass(
  "xcmsScan", 
  slots = c(
    meta = "data.frame",
    intensity = "data.frame",
    filter = "list"
  )
)

setMethod("initialize", "xcmsScan", function(.Object, meta, intensity, filter) {
  
  object <- .Object
  
  if (missing(meta))
    object@meta = data.frame(
      "scanNum" = integer(0), 
      "acquisitionNum" = integer(0), 
      "rt" = numeric(0), 
      "tic" = numeric(0), 
      "peakCount" = integer(0), 
      "msLevel" = integer(0),
      "fromFile" = integer(0), 
      "precScanNum" = integer(0), 
      "precMz" = numeric(0), 
      "precCharge" = integer(0), 
      "precIntensity" = numeric(0), 
      "ID" = character(0),
      stringsAsFactors = FALSE
    ) else
      object@meta <- meta
    
    if (missing(intensity))
      object@intensity = data.frame(
        "mz" = numeric(0), 
        "intensity" = numeric(0), 
        "ID" = character(0),
        stringsAsFactors = FALSE
      ) else
        object@intensity <- intensity
      
      if (missing(filter))
        object@filter <- list() else
          object@filter <- filter
        
        validObject(object)
        return(object)
})


setValidity("xcmsScan", function(object) {
  
  .validRestrictedDataFrame(
    df = object@meta, 
    name = .xcmsViewerInternalObjects$xcmsScan_meta_column,
    class = .xcmsViewerInternalObjects$xcmsScan_meta_column_class,
    str = "object@meta"
  )
  
  .validRestrictedDataFrame(
    df = object@intensity, 
    name = .xcmsViewerInternalObjects$xcmsScan_intensity_column,
    class = .xcmsViewerInternalObjects$xcmsScan_intensity_column_class,
    str = "object@intensity"
  )
  
  if (any(!fastmatch::"%fin%"(object@intensity$ID, object@meta$ID)))
    return("There are scan IDs in xcmsScan@intensity that are not recorded in xcmsScan@meta!")
  
  return(TRUE)
})

################################ xcmsAnnot ############################
setClass(
  "xcmsAnnot", 
  slots = c(
    MS1Annot = "data.frame",
    MS1Meta = "data.frame",
    MS2RefList = "list",
    MS2RefMeta = "data.frame",
    MS2Match = "data.frame"
  ))

# xcmsAnnot
.xcmsViewerInternalObjects$xcmsAnnot_MS1Annot_column <- c(
  "feature_ID", "annot_internal_ID", "adduct", "monoisotopic_weight", "query_weight", 
  "mass_diff", "delta_PPM", "annot_MS2")
.xcmsViewerInternalObjects$xcmsAnnot_MS1Annot_column_class <- c(
  "character", "character", "character", "numeric", "numeric", "numeric", "numeric", "factor")


.xcmsViewerInternalObjects$xcmsAnnot_MS1Meta_column <- c(
  "annot_internal_ID", "name", "chemical_formula", "molecular_weight", 
  "annot_database", "CAS_reg", "InChI", "SMILE", "annot_database_ID")
.xcmsViewerInternalObjects$xcmsAnnot_MS1Meta_column_class <- c(
  "character", "character", "character", "numeric", "character", "character", 
  "character", "character", "character")

.xcmsViewerInternalObjects$xcmsAnnot_MS2Match_column <- c(
  "ms2_scan_ID", "annot_ms2_internal_ID", "cos")
.xcmsViewerInternalObjects$xcmsAnnot_MS2Match_column_class <- c(
  "character", "character", "numeric")

.xcmsViewerInternalObjects$xcmsAnnot_MS2RefMeta_column <- c(
  "annot_ms2_internal_ID", "annot_internal_ID", "collision_energy_level", 
  "collision_energy_voltage", "ionization_mode", "predicted", "database")
.xcmsViewerInternalObjects$xcmsAnnot_MS2RefMeta_column_class <- c(
  "character", "character", "factor", "numeric", "factor", "logical", "character")


setMethod("initialize", "xcmsAnnot", function(
  .Object, MS1Annot, MS1Meta, MS2RefList, MS2RefMeta, MS2Match
) {
  
  # a table for the mapping of ms1
  if (missing(MS1Annot))
    MS1Annot <- data.frame(
      feature_ID = character(0), 
      annot_internal_ID = character(0),
      adduct = character(0),
      monoisotopic_weight = numeric(0),
      query_weight = numeric(0),
      mass_diff = numeric(0),
      delta_PPM = numeric(0),
      annot_MS2 = factor(levels = c("+", "")),
      stringsAsFactors = FALSE
    ) 
  
  # a meta table for ms1
  if (missing(MS1Meta))
    MS1Meta <- data.frame(
      annot_internal_ID = character(0), # PK
      name = character(0),
      chemical_formula = character(0), 
      molecular_weight = numeric(0),
      annot_database = character(0), 
      CAS_reg = character(0),
      InChI = character(0),
      SMILE = character(0),
      annot_database_ID = character(0),
      stringsAsFactors = FALSE
    )
  
  # a table for the mapping between measured ms2 and references ms21
  if (missing(MS2Match))
    MS2Match <- data.frame(
      ms2_scan_ID = character(0),
      annot_ms2_internal_ID = character(0),
      cos = numeric(0),
      stringsAsFactors = FALSE
    )
  
  # a list of reference ms2 spectra
  if (missing(MS2RefList))
    MS2RefList <- list(
      # ms2_reference_1 = data.frame(
      #   mz = numeric(0), 
      #   intensity = numeric(0)
      # )
    )
  
  # a meta table for the reference ms2, linked to the ms1 entries
  if (missing(MS2RefMeta))
    MS2RefMeta <- data.frame(
      annot_ms2_internal_ID = character(0),
      annot_internal_ID = character(0),
      collision_energy_level  = factor(levels = c("low", "med", "high")),
      collision_energy_voltage = numeric(0),
      ionization_mode = factor(levels = c("positive", "negative")),
      predicted = logical(0),
      database = character(0),
      stringsAsFactors = FALSE
    )
  
  .Object@MS1Annot <- MS1Annot
  .Object@MS1Meta <- MS1Meta
  .Object@MS2RefList <- MS2RefList
  .Object@MS2RefMeta <- MS2RefMeta
  .Object@MS2Match <- MS2Match
  
  validObject(.Object)
  return(.Object)
})


setValidity("xcmsAnnot", function(object) {
  
  .validRestrictedDataFrame(
    df = object@MS1Annot,
    name = .xcmsViewerInternalObjects$xcmsAnnot_MS1Annot_column,
    class = .xcmsViewerInternalObjects$xcmsAnnot_MS1Annot_column_class,
    str = "object@MS1Annot"
  )
  .validRestrictedDataFrame(
    df = object@MS2Match,
    name = .xcmsViewerInternalObjects$xcmsAnnot_MS2Match_column,
    class = .xcmsViewerInternalObjects$xcmsAnnot_MS2Match_column_class,
    str = "object@MS2Match"
  )
  .validRestrictedDataFrame(
    df = object@MS2RefMeta,
    name = .xcmsViewerInternalObjects$xcmsAnnot_MS2RefMeta_column,
    class = .xcmsViewerInternalObjects$xcmsAnnot_MS2RefMeta_column_class, 
    str = "object@MS2RefMeta",
    contain = TRUE
  )
  .validRestrictedDataFrame(
    df = object@MS1Meta,
    name = .xcmsViewerInternalObjects$xcmsAnnot_MS1Meta_column,
    class = .xcmsViewerInternalObjects$xcmsAnnot_MS1Meta_column_class, 
    str = "object@MS2RefMeta",
    contain = TRUE
  )
  
  # named list of data.frame
  # the name can be mapped to MS1Annot$annot_database
  # the data.frame has at list one column, named as "annot_database_ID"
  
  
  
})


