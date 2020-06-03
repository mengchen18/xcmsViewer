library(xcmsViewerApp)
library(shinythemes)

dir <- "" 
ff <- list.files(dir, pattern = ".RDS|.rds|.Rds")
vv <- gsub(".RDS|.rds|.Rds", "", ff, ignore.case= TRUE)
names(ff) <- vv

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
  featureSelected <- callModule(
    xcmsViewerApp:::featureStatsTab, id = "fst", dat = dat, dataChanged = reactive(input$selectData) 
    )
  
  #   The left panel
  callModule(
    xcmsViewerApp:::phenoIntensityEIC, id = "pies", dat = dat, featureSelected = featureSelected
    )
  
  #   The right panel - annotation table panel
  callModule(
    xcmsViewerApp:::msAnnotation, id = "msat", dat = dat, featureSelected = featureSelected
    )
  
  #   The right panel - the peak panel
  callModule(
    xcmsViewerApp:::peakEIC, id = "pkeic", dat = dat, featureSelected = featureSelected
    )
  
  #   Free chromatogram
  callModule(xcmsViewerApp:::freeChrom, id = "fc", dat = dat)
  
  session$onSessionEnded(function() {
    stopApp()
  })
}
