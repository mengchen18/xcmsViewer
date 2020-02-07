freeChromUI <- function(id) {
  ns <- NS(id)  
  
  sidebarLayout(
    sidebarPanel(
      tags$div(
        style="margin-bottom:8px;",
        DT::dataTableOutput(ns("pheno"))
      ),
      fluidRow(
        column(
          6,
          checkboxInput(ns("multiSelect"), "Mutiple selection", value = TRUE)
        ),
        column(
          6, align = 'right',
          actionButton(ns("clearRow"), "Clear selection")
        )
      )
    ),
    mainPanel(
      wellPanel(
        fluidRow(
          column(
            3, textInput(ns("range_RT"), label = "Range RT", placeholder = "e.g. 120-130")
          ),
          column(
            3, textInput(ns("range_MZ"), label = "Range MZ", placeholder = "e.g. 320-321 or 455")
          ),
          column(
            4, selectizeInput(ns("phenoCols"), "Color:", choices = NULL)
          ),
          column(
            2, 
            tags$div(
              style="margin-top:25px;",
              actionButton(ns("sbutton"), "Apply")
            )
          )
        ),
        plotEIC_UI(ns("eic_plot"), height = 325)
      )
    )
  )
}


freeChrom <- function(input, output, session, dat) {
  
  output$pheno <- DT::renderDataTable({
    
    if (input$multiSelect)
      ss <- "multiple" else
        ss <- "single"
      
      DT::datatable(
        dat()$pheno,
        selection = ss,
        filter = "top",
        rownames = FALSE,
        option = list(scrollX = TRUE, scrollY = "256px", dom = "t")
      )
  })
  
  dtProxy <- DT::dataTableProxy('pheno')
  observeEvent(input$clearRow, {
    DT::selectRows(dtProxy, NULL)
  })
  
  observe(
    updateSelectizeInput(session, "phenoCols", choices = colnames(dat()$pheno))
    )
  
  range <- eventReactive(input$sbutton, {
    
    if (input$range_RT == "" & input$range_MZ == "") {
      return(NULL)
    } else if (input$range_RT != "" & input$range_MZ == "") {
      l <- list( 
        rt = parseRange(input$range_RT),
        mz = c(-Inf, Inf) 
      )
    } else if (input$range_RT == "" & input$range_MZ != "") {
      l <- list( 
        rt = c(-Inf, Inf),
        mz = parseRange(input$range_MZ)
      ) 
    } else {
      l <- list(
        rt = parseRange(input$range_RT),
        mz = parseRange(input$range_MZ)
      ) 
    }
    l
  })
  
  x <- reactive({
    req(range())
    r <- eic(
      itab = dat()$scanIntensityTab, 
      mtab = dat()$scanMetaTab, 
      rt = range()$rt, 
      mz = range()$mz
    )
    req(r)
    r
  })
  
  colvar <- callModule(
    plotEIC_module, "eic_plot",
    react_x = x,
    react_diffPheno = reactive( 
      as.character(dat()$pheno[[input$phenoCols]]) 
    ),
    react_select = reactive( 
      input$pheno_rows_selected 
    )
  )
}

# library(shiny)
# 
# ui <- fluidPage(
#   freeChromUI("fc")
# )
# server <- function(input, output, session) {
#   callModule(freeChrom, id = "fc")
# }
# 
# shinyApp(ui, server)

