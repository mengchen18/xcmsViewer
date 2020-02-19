plotEIC_UI <- function(id, height = 300) {
  ns <- NS(id)
  plotOutput(
    ns("plot1"), height = height, dblclick = ns("plot1_dblclick"),
    brush = brushOpts(
      id = ns("plot1_brush"), resetOnNew = TRUE
    ),
    click = ns("eic_click")
  )
}

plotEIC_module <- function(
  input, output, session, 
  react_x, react_diffPheno = reactive(NULL), 
  react_select = reactive(NULL), react_vline = reactive(NULL),
  react_initXlim = reactive(NULL), 
  reset_click = reactive(NULL)
) {

  colorLegend <- reactive({
    req(react_x())
    req(react_diffPheno())
    colCode(x = react_x(), diffPheno=react_diffPheno(), select = react_select())
    })
  rtrange <- reactive(
    range(react_x()$rt)
    )
  
  ranges <- reactiveValues(x=NULL, y=NULL)
  observe({    
    if (!is.null(react_initXlim())) {
      ranges$x <- react_initXlim()
      ranges$y <- NULL    
    }
  }) 
  
  output$plot1 <- renderPlot({
    plotEIC(x=react_x(), unifile = colorLegend()$unifile, col = colorLegend()$col,
            diffPheno = colorLegend()$diffPheno, vline = react_vline(), 
            xlim = ranges$x, ylim = ranges$y)    
  })

  observeEvent(react_x(), {
    ranges$x <- NULL
    ranges$y <- NULL
    if (!is.null(react_initXlim()))
      ranges$x <- react_initXlim()
  })
  
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
      if (!is.null(react_initXlim()))
        ranges$x <- react_initXlim()
    }
  })

  file_click <- reactiveVal(NULL)
  observeEvent(input$eic_click, {
    req(input$eic_click)
    np <- nearPoints(
      react_x(), input$eic_click, xvar = "rt", yvar = "intensity", maxpoints = 1, threshold=10
      )
    file_click(np$file)
    })
  observeEvent(reset_click(), {
    file_click(NULL)
    })

  reactive(
    c(colorLegend(), list(file_clicked = file_click()))
    )
  }


# ui <- fluidPage(
#   sidebarLayout(
#     sidebarPanel(),
#     mainPanel(
#       plotEIC_UI("test")
#     )
#   )
# )
# 
# server <- function(input, output, session) {
#   callModule(plotEIC_module, "test", react_x = reactive(x), react_diffPheno = reactive(c("T1", "T1", "T3", "T3")))
#   # callModule(plotEIC_module, "test", react_x = reactive(qq))
# }
# #
# shinyApp(ui, server)

