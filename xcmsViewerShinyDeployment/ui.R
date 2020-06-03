library(xcmsViewerApp)
library(shinythemes)

## change here
dir <- "" 
ff <- list.files(dir, pattern = ".RDS|.rds|.Rds")
vv <- gsub(".RDS|.rds|.Rds", "", ff, ignore.case= TRUE)
names(ff) <- vv

ui = navbarPage(
  "XCMSViewer - Explore Untargeted Metaboloimcs Data",
  theme = shinytheme("yeti"),
  
  tabPanel(
    "Metaboloimc feature",
    fluidRow(
      # left panel - phenotype related
      column( 3, xcmsViewerApp:::phenoIntensityEIC_UI("pies") ),
      # mid panel - table
      column( 6, xcmsViewerApp:::featureStatsTab_UI("fst") ),
      # right panel - annotation relalted
      column( 3,
              tabsetPanel(
                tabPanel( "Annotation", xcmsViewerApp:::msAnnotation_UI("msat") ),
                tabPanel( "Peaks", xcmsViewerApp:::peakEIC_UI("pkeic") )
              )
      )
    )
  ),
  tabPanel(
    "Free chromatogram",
    xcmsViewerApp:::freeChromUI("fc")
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
)

