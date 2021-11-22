#' Start xcmsViewer app
#' @param dir directory contain processed data as .RDS format, can be load into xcmsViewer
#' @param ... current not used. In previous version, it was for other arguments passed to shiny::runApp
#' @export
#' @rawNamespace import(shiny, except = span)
#' @importFrom shinythemes shinytheme
#' @importFrom beeswarm beeswarm
#' @importFrom grDevices rgb
#' @importFrom utils packageVersion
#' @import RSQLite DBI
#' @import ExpressionSetViewer

xcmsViewer <- function(dir, ...) {
  otherParams <- list(...)
  
  ll <- list(
    list(
      tabName = "Annotation",
      moduleName = "xcmsannot",
      moduleUi = xcmsAnnotationTab_ui,
      moduleServer  = xcmsAnnotationTab_module
    )
    # # This function is disabled for the current version
    # # perhaps will be restored in the future
    # ,
    # list(
    #   tabName = "Free Chromatogram",
    #   moduleName = "fChrom",
    #   moduleUi = freeChromUI,
    #   moduleServer  = freeChrom
    # )
  )
  returnValueOrNULL <- function(x, v) {
    r <- x$header[x$axis == v]
    if (is.na(r))
      return(NULL)
    r
  }

  ExpressionSetViewer::ExpressionSetViewer(
    dir,
    filePattern = ".RDS$|.DB$",
    esetLoader = function(x) {
      if (grepl(".db$", x, ignore.case = TRUE)) {        
        obj <- dbConnect(RSQLite::SQLite(), x)
        x <- RSQLite::dbGetQuery(obj, 'SELECT * FROM defaultVis;')
        attr(obj, "fx") <- returnValueOrNULL(x, "fx")
        attr(obj, "fy") <- returnValueOrNULL(x, "fy")
        attr(obj, "sx") <- returnValueOrNULL(x, "sx")
        attr(obj, "sy") <- returnValueOrNULL(x, "sy")
      } else if (grepl("RDS", x, ignore.case = TRUE)) {
        obj <- readRDS(x)
      }
      obj
      },
    exprsGetter = function(x)  {
      getFeatureExprs(x)      
    },
    pDataGetter = function(x) {
      getPheno(x)
    },
    fDataGetter = function(x) {
      getFeatureMeta(x)
    },
    additionalTabs = ll,
    appName = "XcmsViewer",
    appVersion = packageVersion("xcmsViewer")
  )
  
}

