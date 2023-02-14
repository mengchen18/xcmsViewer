
landingPage_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style("#landingpageid{margin: auto;}"),
    h1("Welcome to BayBioMS xcmsViewer for untargeted metabolomics data!"),
    absolutePanel(id = "landingpageid", class = "myClass", fixed = TRUE
                  ,draggable = FALSE, left = 0, right = 0, top = "35%"
                  ,width = 500, height = "auto", 
                  wellPanel(
                    uiOutput(ns("error.ui")),
                    fluidRow(
                      column(10, passwordInput(ns("id"), "Your passcode", placeholder = "Enter your passcode here!")),
                      column(2, actionButton(ns("go"), "Go!"), style = "padding-top:25px")
                    )
                  )
    )
  )
}

landingPage_module <- function(input, output, session) {
  
  ns <- session$ns
  
  hashTab <- read.delim(
    # "/dss/ga52cib/xcmsViewer/passcodeList.tsv", header = TRUE,
    "/home/shiny/data/xcmsViewer/passcodeList.tsv", header = TRUE, stringsAsFactors = FALSE
    # "/media/share_LRZCC/xcmsViewer/passcodeList_debug.tsv", header = TRUE, stringsAsFactors = FALSE
  )
  
  # hashTab <- data.frame(
  #   code = 666,
  #   path = "~/xcmsViewerDBs/",
  #   stringsAsFactors = FALSE
  # )
  
  m0 <- eventReactive(input$go, {
    id <- trimws(input$id)
    i <- which(hashTab[, 1] %in% id)
    if (length(i) != 1)
      return(NULL)
    i
  })
  
  output$error <- renderUI({
    HTML("There was an error with your passcode, please try again! <br><hr>")
  })
  
  output$error.ui <- renderUI({
    req(is.null(m0()))
    uiOutput(ns("error"))
  })
  
  reactive({
    req(m0())
    hashTab[m0(), 2]
  })
}
