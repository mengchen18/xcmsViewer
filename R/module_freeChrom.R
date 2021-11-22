freeChromUI <- function(id) {
  ns <- NS(id)  
  EIC_ui(id = ns("eic_plot"), layoutFun = function(EICplot, colorControl, repPeakControl) {
    ###
    tagList(
      fluidRow(
        column(
          4, textInput(ns("range_RT"), label = "Range RT", placeholder = "e.g. 120-130")
        ),
        column(
          4, textInput(ns("range_MZ"), label = "Range MZ", placeholder = "e.g. 320-321 or 455")
        ),
        column(
          1, 
          tags$div(
            style="margin-top:25px;",
            actionButton(ns("sbutton"), "Apply")
          )
        ),
        column(
          3, align = "right",
          tags$div(
            style="margin-top:25px;",
            shinyWidgets::dropdownButton(
              circle = FALSE, label = "highlight sample", right = TRUE, width = 800,
              colorControl,
              DT::dataTableOutput(ns("pheno")),
              fluidRow(
                column(
                  6, checkboxInput(ns("multiSelect"), "Mutiple selection", value = TRUE)
                ),
                column(
                  6, actionButton(ns("clearRow"), "Clear selection"), align = "right"
                )
              )
            )
          )
        ),
        column(
          12,
          EICplot,
          conditionalPanel("1 == 2", repPeakControl)
        )
      )
    )
  })
}


freeChrom <- function(
  input, output, session, pdata, fdata, expr, feature_selected, sample_selected, object
) {
  
  output$pheno <- DT::renderDataTable({
    
    if (input$multiSelect)
      ss <- "multiple" else
        ss <- "single"
      
      DT::datatable(
        pdata(),
        selection = ss,
        filter = "top",
        rownames = FALSE,
        option = list(scrollX = TRUE, scrollY = "225px", dom = "t")
      )
  })
  
  dtProxy <- DT::dataTableProxy('pheno')
  observeEvent(input$clearRow, {
    DT::selectRows(dtProxy, NULL)
  })
  
  observe(
    updateSelectizeInput(session, "phenoCols", choices = colnames(pdata()))
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
    list(eic_param = l)
  })
  
  hl <- reactiveVal(NULL)
  observe({
    if (length(s <- sample_selected()) > 0)
      hl(match(s, rownames(pdata())))
  })
  observe({
    if (length(s <- input$pheno_rows_selected) > 0)
      hl(s)
  })
  callModule(EIC_module, id = "eic_plot", pdata = pdata, object = object, obj = range, switchOnRepPeak = Inf, 
             file = hl)
}
