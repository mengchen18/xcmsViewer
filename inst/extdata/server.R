library(xcmsViewer)
library(omicsViewer)
source("/home/shiny/app/landingPage.R")
# source("landingPage.R")

ll <- list(
  list(
    tabName = "Annotation",
    moduleName = "xcmsannot",
    moduleUi = xcmsViewer:::xcmsAnnotationTab_ui,
    moduleServer  = xcmsViewer:::xcmsAnnotationTab_module
  )
)

returnValueOrNULL <- function(x, v) {
  r <- x$header[x$axis == v]
  if (is.na(r))
    return(NULL)
  r
}

server <- function(input, output, session) {
  
  ns <- session$ns
  
  v <- callModule(landingPage_module, id = "test")
  showLanding <- reactiveVal(TRUE)
  
  observe({
    showLanding ( length(v()) != 1 )
  })
  
  output$uis <- renderUI({
    req(showLanding())
    landingPage_ui("test") 
  })
  
  observe(print(v()))
  
  callModule(
    omicsViewer:::app_module, id = "app", .dir = v,
    filePattern = ".RDS$|.DB$",
    esetLoader = function(x) {
      if (grepl(".db$", x, ignore.case = TRUE)) {        
        obj <- RSQLite::dbConnect(RSQLite::SQLite(), x)
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
    exprsGetter = getFeatureExprs,
    pDataGetter = getPheno,
    fDataGetter = getFeatureMeta,
    defaultAxisGetter = function(x, what) attr(x, what),
    additionalTabs = ll,
    appName = "XcmsViewer",
    appVersion = packageVersion("xcmsViewer")
  )

  output$aout <- renderUI({
    req(!showLanding())
    omicsViewer:::app_ui("app")
  })
}
