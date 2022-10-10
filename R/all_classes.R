############################# global variables ########################
.xcmsViewerInternalObjects <- function() {
  list(
    # xcmsFeatureSet
    xcmsFeatureSet_fdata_name = c(
      "mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", "npeaks", 
      "peakidx", "ID", "annot_ms1", "annot_ms2", "ms2spectrum", "purity"),
    xcmsFeatureSet_fdata_name_class = c(
      "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
      "numeric", "list", "character", "character", "character", "character", "numeric"),
    
    # xcmsScan
    xcmsScan_meta_column = c(
      "scanNum", "acquisitionNum", "rt", "tic", "peakCount", "msLevel", "fromFile", 
      "precScanNum", "precMz", "precCharge", "precIntensity", "ID", "peakCountFiltered", "validMS2"),
    xcmsScan_meta_column_class = c(
      "integer", "integer", "numeric", "numeric", "integer", "integer", "integer", 
      "integer", "numeric", "integer", "numeric", "character", "integer", "logical"),
    xcmsScan_intensity_column = c("mz", "intensity", "ID"),
    xcmsScan_intensity_column_class = c("numeric", "numeric", "character"),
    
    # xcmsPeak
    xcmsPeak_table_column = c(
      "ID", "sample", "masterPeak", "into", "mz", "mzmin", 
      "mzmax", "rt", "rtmin", "rtmax", "sn", "is_filled", "ms_level", "ms2Scan", "validMS2" ),
    xcmsPeak_table_column_class = c(
      "character", "integer", "character", "numeric", "numeric", "numeric", 
      "numeric", "numeric", "numeric", "numeric", "numeric", "logical", "integer", "character", "logical"),
    
    # xcmsAnnot
    xcmsAnnot_column = c("ID", "InChIKey", "CID", "cpdName", "formula", "monoMass", "Adduct", 
      "MassDiff", "Nmol", "IPS", "MassWithAdduct", "MassQueried", "DeltaPPM", "smiles", "RT", 
      "ms2_mass", "ms2_intensity", "ms2_purity", "ms2_cos", "score_ms1", "score_ms2", 
      "score_rt", "Score"),
    xcmsAnnot_column_class =  c("character", "character", "character", "character", "character", 
      "numeric", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
      "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric", 
      "numeric", "numeric")
  )
}

############################ new classes ############################

#' @import fastmatch
#' @import Biobase
#' @import xcms
#' @import methods
#' @importClassesFrom xcms PeakGroupsParam
#' @importClassesFrom xcms ObiwarpParam
#' @importClassesFrom xcms PeakDensityParam
#' @importClassesFrom xcms MzClustParam
#' @importClassesFrom xcms NearestPeaksParam
#' @importClassesFrom xcms CentWaveParam
#' @importClassesFrom xcms CentWavePredIsoParam
#' @importClassesFrom xcms MatchedFilterParam
#' @importClassesFrom xcms MassifquantParam
#' @importClassesFrom xcms MSWParam

setClassUnion(
  name = "xcmsRTAdjustParamClass", 
  members = c("PeakGroupsParam", "ObiwarpParam", "NULL")
)

setClassUnion(
  name = "xcmsPeakGroupParamClass",
  c("PeakDensityParam", "MzClustParam", "NearestPeaksParam", "NULL")
)

setClassUnion(
  name = "listORnull",
  c("list", "NULL")
)

setClassUnion(
  "chromPeakParamORnull", 
  c("NULL",
    "CentWaveParam",
    "CentWavePredIsoParam",
    "MatchedFilterParam",
    "MassifquantParam",
    "MSWParam"))

setClassUnion("numericORnull", c("NULL", "numeric"))

############################ xcmsFeatureSet ############################
#' xcmsFeatureSet S4 class that extends \code{"ExpressionSet"} class.
#'
#' Some details 
#'
#' @slot RTAdjustParam Parameter used for retention time alignment
#' @slot peakGroupParam Parameter used for peak grouping
#'
#' @name xcmsFeatureSet-class
#' @rdname xcmsFeatureSet-class
#' @export
#' 
setClass(
  "xcmsFeatureSet", 
  slots = c(
    "RTAdjustParam" = "xcmsRTAdjustParamClass", 
    "peakGroupParam" = "xcmsPeakGroupParamClass"
  ),
  contains = "ExpressionSet"
)

#' initialize xcmsFeatureSet
#' description 3
#' @param .Object object 
#' @param intensity intensity 
#' @param masterPeak maste peak
#' @param RTAdjustParam rt adjustment time
#' @param peakGroupParam peak group
#' @param ... other parameter
#' @rdname xcmsFeatureSet-class
#' 
setMethod("initialize", "xcmsFeatureSet", function(
  .Object, intensity, masterPeak, RTAdjustParam=NULL, peakGroupParam=NULL,
  ...) {
  object <- .Object
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
  obj@RTAdjustParam <- RTAdjustParam
  obj@peakGroupParam <- peakGroupParam
  validObject(obj)
  obj
})

setValidity("xcmsFeatureSet", function(object) {
  pd <- fData(object)
  .validRestrictedDataFrame(
    pd,
    name = .xcmsViewerInternalObjects()$xcmsFeatureSet_fdata_name,
    class = .xcmsViewerInternalObjects()$xcmsFeatureSet_fdata_name_class,
    str = "xcmsFeatureSet@fData",
    contain = TRUE
  )
})



############################ xcmsPeak ############################
#' xcmsPeak S4 class 
#'
#' Peak tables extracted using xcms peakpicking
#'
#' @slot table peak table
#' @slot param Parameter used for peak picking
#' @slot postFilter Parameter used for post filtering of peaks
#' 
#' @name xcmsPeak-class
#' @rdname xcmsPeak-class
#' @export
#' 
setClass(
  "xcmsPeak", 
  slots = c(
    table = "data.frame",
    param = "chromPeakParamORnull",
    postFilter = "numericORnull"
  )
)

#' initialize xcmsPeak
#' description 2
#' @param .Object object
#' @param table table
#' @param param param
#' @param postFilter postFilter
#' @rdname xcmsPeak-class
#' 
setMethod("initialize", "xcmsPeak", function(.Object, table, param, postFilter) {
  object <- .Object
  if (missing(table)) {
    object@table = data.frame(
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
    object@table <- table[.xcmsViewerInternalObjects()$xcmsPeak_table_column]
  
  if (missing(param))
    object@param <- NULL else
      object@param <- param
    
    if (missing(postFilter))
      object@postFilter <- NULL else
        object@postFilter <- postFilter
      
      validObject(object)
      return(object)
      
})

setValidity("xcmsPeak", function(object) {
  .validRestrictedDataFrame(
    df = object@table, 
    name = .xcmsViewerInternalObjects()$xcmsPeak_table_column, 
    class = .xcmsViewerInternalObjects()$xcmsPeak_table_column_class,
    str = "object@table"
  )
})


############################ xcmsScan ############################
#' xcmsScan S4 class 
#'
#' scans
#'
#' @slot meta scan meta
#' @slot intensity The mz and intensities of scans
#' @slot filter the filter used to filter the scans
#' @slot keepMS1 logical; whether MS1 scans are kept
#' 
#' @name xcmsScan-class
#' @rdname xcmsScan-class
#' @export
#' 
setClass(
  "xcmsScan", 
  slots = c(
    meta = "data.frame",
    intensity = "data.frame",
    filter = "list",
    keepMS1 = "logical"
  )
)

#' initialize xcmsScan
#' description 1
#' @param .Object object
#' @param meta meta
#' @param intensity intensity
#' @param filter filter
#' @param keepMS1 logical; whether MS1 scans should be kept
#' @rdname xcmsScan-class
#' 
setMethod("initialize", "xcmsScan", function(.Object, meta, intensity, filter, keepMS1 = TRUE) {
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
        
        object@keepMS1 <- keepMS1
        validObject(object)
        return(object)
})

setValidity("xcmsScan", function(object) {
  
  .validRestrictedDataFrame(
    df = object@meta, 
    name = .xcmsViewerInternalObjects()$xcmsScan_meta_column,
    class = .xcmsViewerInternalObjects()$xcmsScan_meta_column_class,
    str = "object@meta"
  )
  
  .validRestrictedDataFrame(
    df = object@intensity, 
    name = .xcmsViewerInternalObjects()$xcmsScan_intensity_column,
    class = .xcmsViewerInternalObjects()$xcmsScan_intensity_column_class,
    str = "object@intensity"
  )
  
  if (any(!fastmatch::"%fin%"(object@intensity$ID, object@meta$ID)))
    return("There are scan IDs in xcmsScan@intensity that are not recorded in xcmsScan@meta!")
  
  return(TRUE)
})


################################# prunedXcmsSet ################################
#' prunedXcmsSet S4 class 
#'
#' prunedXcmsSet
#'
#' @slot featureSet xcmsFeatureSet
#' @slot peak xcmsPeak
#' @slot scan xcmsScan
#' @slot annot annotation
#' @slot EIC extracted ion currents
#' 
#' @name prunedXcmsSet-class
#' @rdname prunedXcmsSet-class
#' @export
#' 
setClass(
  "prunedXcmsSet", 
  slot = c(
    featureSet = "xcmsFeatureSet",
    peak = "xcmsPeak",
    scan = "xcmsScan", 
    annot = "data.frame",
    EIC = "list"
  )
)

setValidity(Class = "prunedXcmsSet", function(object) {
  validObject(object@featureSet)
  validObject(object@peak)
  validObject(object@scan)
  .validRestrictedDataFrame(
    df = object@annot, 
    name = .xcmsViewerInternalObjects()$xcmsAnnot_column,
    class = .xcmsViewerInternalObjects()$xcmsAnnot_column_class,
    str = "object@intensity"
  )
})

