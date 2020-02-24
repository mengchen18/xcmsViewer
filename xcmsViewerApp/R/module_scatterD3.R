scatterD3_UI <- function(id) {
  ns <- NS(id)
  tagList(
    div(style="display: inline-block;vertical-align:top; width: 24%;", varSelectInput(ns("sctd3_var_x"), label = "X-axis", data = NULL) ),
    div(style="display: inline-block;vertical-align:top; width: 24%;", varSelectInput(ns("sctd3_var_y"), label = "Y-axis", data = NULL) ),    
    div(style="display: inline-block;vertical-align:top; width: 24%;", varSelectInput(ns("sctd3_var_col"), label = "Point color", data = NULL) ),   
    div(style="display: inline-block;vertical-align:top; width: 24%;", sliderInput(ns("sctd3_var_p_opacity"), label = "Point opacity", min = 0, max = 1, value = 0.7)), 
    div(style="display: inline-block;vertical-align:top; width: 24%;", varSelectInput(ns("sctd3_var_pch"), label = "Point shape", data = NULL) ),    
    div(style="display: inline-block;vertical-align:top; width: 24%;", varSelectInput(ns("sctd3_var_cex"), label = "Point size", data = NULL) ),    
    div(style="display: inline-block;vertical-align:top; width: 24%;", varSelectInput(ns("sctd3_var_tooltips"), label = "Tooptips Info", data = NULL, multiple = TRUE) ),
    div(style="display: inline-block;vertical-align:top; width: 24%;", sliderInput(ns("sctd3_var_sizerange"), label = "Point size range", min = 1, max = 1000, value = c(20, 400)) )
  )
}

scatterD3_plot <- function(id) { 
  ns <- NS(id)
  scatterD3::scatterD3Output(ns("scatterPlot"))
}

scatterD3_Module <- function(  
  input, output, session, data, rows=reactive(NULL), clickRetVar = reactive(NULL) 
) {

  data2 <- reactive({
    if (is.null(rows()))
      return(data())
    data()[rows(), ]
    })
  
  f0 <- function(x, v) {
    if (length(xx <- intersect(v, x)) > 0)
      return(xx)
    NA
  }

  observe( updateVarSelectizeInput(session, "sctd3_var_x", data = data(), selected = f0("rtmed", colnames(data()))) )
  observe( updateVarSelectizeInput(session, "sctd3_var_y", data = data(), selected = f0("mzmed", colnames(data())) ) )
  observe( updateVarSelectizeInput(session, "sctd3_var_col", data = data(), selected = f0("QC", colnames(data()))) )
  observe( updateVarSelectizeInput(session, "sctd3_var_pch", data = data(), selected = NA) )
  observe( updateVarSelectizeInput(session, "sctd3_var_cex", data = data(), selected = f0("npeaks", colnames(data()))) )
  observe( updateVarSelectizeInput(session, "sctd3_var_tooltips", data = data(), selected = NA) )
  
  output$scatterPlot <- scatterD3::renderScatterD3({

    req(data2())
    req(input$sctd3_var_x)
    req(input$sctd3_var_y)
    fn <- function(x) {
      if (length(x) == 0)
        return(NULL)
      x }
    
    var_x <- as.character(input$sctd3_var_x)
    var_y <- as.character(input$sctd3_var_y)
    var_col <- as.character(input$sctd3_var_col)
    var_pch <- as.character(input$sctd3_var_pch)
    var_cex <- as.character(input$sctd3_var_cex)

    tv <- as.character(input$sctd3_var_tooltips)
    if (isTruthy(tv)) {
      ttvar <- lapply( tv, function(x) paste(sprintf("<strong>%s</strong>", x), data2()[, x]) )
      ttvar <- do.call(function(...) paste(..., sep="<br />"), ttvar)
      print(ttvar)
    } else
      ttvar <- NULL    

    scatterD3::scatterD3(
      x = data2()[, var_x], y = data2()[, var_y],
      
      transitions =TRUE,
      key_var = fn( data2()[, clickRetVar()] ), 
      size_range = input$sctd3_var_sizerange, 
      hover_size = 1.5,
      hover_opacity = 1,
      tooltip_text = ttvar,
      point_opacity = input$sctd3_var_p_opacity, 
      col_var = fn( data2()[, var_col] ),
      col_lab = var_col,
      symbol_var = fn( data2()[, var_pch] ),
      symbol_lab = var_pch,
      size_var = fn( data2()[, var_cex] ),
      size_lab = var_cex,
      xlab = var_x, 
      ylab = var_y,
      click_callback = sprintf("function(id, index) {
              if(id && typeof(Shiny) != 'undefined') {
                Shiny.onInputChange('%s', index);}
              }", session$ns("selected_point"))
      )
  })

  reactive({
    if (is.null(clickRetVar()))
      return(input$"selected_point") else
        data2()[input$"selected_point", clickRetVar()]
  })
}


# library(shiny)
# library(scatterD3)

# source("/media/share_baybioms/Projects/008_Bioinformatics/B012_MAIT_devel/xcmsViewerApp/R/module_scatterD3.R")
# x  <- mtcars
# x$nam <- make.names(rownames(x))

# ui <- fluidPage(
#   scatterD3_UI("test"),
#   scatterD3_plot("test")
# )

# server <- function(input, output, session) {
#   v <- callModule(
#     scatterD3_Module, id = "test", data = reactive(x), clickRetVar = reactive("nam")
#   )
  
#   observe( print(v()))
# }

# shinyApp(ui, server)