#' Start xcmsViewer app
#' @param dir directory contain processed data as .RDS format, can be load into xcmsViewer
#' @param ... other arguments passed to shiny::runApp
#' @export
#' @import shiny
#' @importFrom shinythemes shinytheme
#' @importFrom beeswarm beeswarm
#' @importFrom grDevices rgb
#' @examples 
#' # library(xcmsViewerApp)
#' # f <- system.file(package = "xcmsViewerApp", "extdata")
#' # xcmsViewer(f)

xcmsViewer <- function(dir, ...) {
  
  ff <- list.files(dir, pattern = ".RDS")
  vv <- gsub(".RDS", "", ff, ignore.case= TRUE)
  names(ff) <- vv
  
  app <- list(
    ui = navbarPage(
      "XCMSViewer - Explore Untargeted Metaboloimcs Data",
      theme = shinytheme("yeti"),
      
      tabPanel(
        "Metaboloimc feature",
        fluidRow(
          # left panel - phenotype related
          column( 3, phenoIntensityEIC_UI("pies") ),
          # mid panel - table
          column( 6, featureStatsTab_UI("fst") ),
          # right panel - annotation relalted
          column( 3,
                  tabsetPanel(
                    tabPanel( "Annotation", msAnnotation_UI("msat") ),
                    tabPanel( "Peaks", peakEIC_UI("pkeic") )
                  )
          )
        )
      ),
      tabPanel(
        "Free chromatogram",
        freeChromUI("fc")
      ),
      absolutePanel(
        style="z-index:1000;",
        selectizeInput("selectData", NULL, choices = vv, 
                       options = list(placeholder = 'Select a dataset',
                                      onInitialize = I('function() { this.setValue(""); }')
                       )
        ),
        right = "20px", top="5px"
      )
    ),
    
    # Define server logic required to draw a histogram
    server = function(input, output, session) {
      
      observe(
        print(sprintf("Data selected: \n%s", input$selectData))
      )
      
      dat <- reactive({        
        req(input$selectData)
        showModal(
          modalDialog(
            title = "Loading data ...",
            tags$h3(input$selectData),
            footer=tags$h6("This may take a few minutes!")
            )
          )
        d <- readRDS(file.path(dir, ff[input$selectData]))
        removeModal()
        d
      })
      
      #   The mid panel - feature table
      featureSelected <- callModule(featureStatsTab, id = "fst", dat = dat, dataChanged = reactive(input$selectData) )
      
      #   The left panel
      callModule(phenoIntensityEIC, id = "pies", dat = dat, featureSelected = featureSelected)
      
      #   The right panel - annotation table panel
      callModule(msAnnotation, id = "msat", dat = dat, featureSelected = featureSelected)
      
      #   The right panel - the peak panel
      callModule(peakEIC, id = "pkeic", dat = dat, featureSelected = featureSelected)
      
      #   Free chromatogram
      callModule(freeChrom, id = "fc", dat = dat)
      
      session$onSessionEnded(function() {
        stopApp()
      })
    }
  )
  
  runApp(app, ...)
}

