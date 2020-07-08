
.validRestrictedDataFrame <- function(df, name, class, str, contain = FALSE) {
  if (contain) {
    i <- setdiff(name, colnames(df))
    if (length(i) > 0) 
      return(sprintf("Missing required columns in %s: %s", str, paste(i, collapse = ", ")))
    df <- df[, name]
  } else {
    if (!identical(colnames(df), name))
      return(sprintf("Problem in %s column name!", str))
  }
  if (!all( sapply(df, class) == class))
    return(sprintf("Problem in %s column class!", str))
}


setGeneric("masterPeak", function(.Object) {
  standardGeneric("masterPeak")
})
setMethod("masterPeak", "xcmsFeatureSet", function(.Object) {
  assayData(.Object)$masterPeak
})
