mirrorPlotUI <- function(id) {
  ns <- NS(id)
  plotOutput(ns("plot1"), click = ns("plot_click"), 
             brush = brushOpts(ns("plot_brush"), resetOnNew = TRUE), 
             dblclick = ns("plot_dblclick"))
}

mirrorPlotModule <- function(input, output, session, measured, standard, ppmtol, legend.measured, legend.standard) {
  
  v <- reactive({
    req(!is.null(measured()) || !is.null(standard()))
    .prep_mirrorPlot(peak.upper=measured(), peak.lower=standard(), ppmtol=ppmtol()) 
  })

  ranges <- reactiveValues(x = NULL, y = NULL)
  
  observeEvent(v(), {
    ranges$x <- NULL
  })
  
  observeEvent(input$plot_dblclick, {
    brush <- input$plot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
    } else {
      ranges$x <- NULL
    }
  })
  
  output$plot1 <- renderPlot({
    .mirrorPlot(v()$tab, v()$col,legend.peak.upper=legend.measured(), 
                legend.peak.lower=legend.standard(), xlim = ranges$x)
  })
  
  output$info <- renderText({
    paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
  })
}

# 
# # ### example
# dat <- readRDS("/media/share_baybioms/Projects/008_Bioinformatics/B012_MAIT_devel/Dat/00_exampleData_processed.RDS")
# 
# set.seed(1000)
# rt <- dat()$scanIntensityTab[dat()$scanIntensityTab$ID %in% "F1.S2345", ]
# me <- rt[sample(1:rnow(rt), size = 5), ] # + rnorm(5, sd = 0.0005)
# st <- rt[sample(1:nrow(rt), size = 8), ]
# 
# ui <- basicPage(
#   mirrorPlotUI("try")
# )
# 
# server <- function(input, output, session) {
#   callModule(mirrorPlotModule, "try", measured = reactive(me),
#              standard = reactive(st), legend.measured  = reactive("xxxxxxxxxxx"),
#              legend.standard = reactive("yyyyyyyyyyyyy"), ppmtol = reactive(10))
# }
# 
# shinyApp(ui, server)
# 
# #####
# rt <- dat()$scanIntensityTab[dat()$scanIntensityTab$ID %in% "F1.S2345", ]
# me <- rt
# st <- me
# 
# st2 <- st
# st2$mz <- st2$mz+4
# mirrorPlot(peak.upper = me, peak.lower = st2)
# mirrorPlot(peak.upper = me, peak.lower = NULL)
# mirrorPlot(peak.upper = me, peak.lower = NULL, xlim = c(50, 60))
