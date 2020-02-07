peakEIC_UI <- function(id) {
	ns <- NS(id)
	shiny::tagList(
		wellPanel(
              DT::dataTableOutput(outputId = ns("peakTable")),
              style = "padding: 5px 15px 5px 15px; font-size: 12px"
            ),
            wellPanel(
              plotEIC_UI(ns("peak_eic"))
            )
		)
}

peakEIC <- function(input, output, session, dat, featureSelected = reactive(NULL)) {
	output$peakTable <- DT::renderDataTable({
    req(featureSelected()$peakTab)
    DT::datatable(
      featureSelected()$peakTab,
      selection = "single",
      rownames = FALSE,
      options = list(scrollX = TRUE)
    )
    })
  # when select a peak, returns
  peakSelected <- reactive({ 
    req(input$peakTable_rows_selected)
    f <- featureSelected()$peakTab[input$peakTable_rows_selected, ]
        
    mzoff <- 25*1e-6*f$mz
    res <- list()
    
    # 1. ms2 scans
    res$ms2scans <- f$ms2Scan  
    
    # 2. peak EIC
    tabrtmin <- min(featureSelected()$peakTab$rtmin)-20
    tabrtmax <- max(featureSelected()$peakTab$rtmax)+20
    x <- eic(
      itab = dat()$scanIntensityTab, 
      mtab = dat()$scanMetaTab,
      rt = c(tabrtmin, tabrtmax), 
      mz = c(f$mzmin - mzoff, f$mzmax + mzoff),
      file = f$sample)
    res$eic <- x
    
    # 3. x range
    res$rtrange <- c(tabrtmin, tabrtmax)
    
    # 4. react_vline
    res$rtVline <- c(
      rtmed = f$rt,
      rtmin = f$rtmin,
      rtmax = f$rtmax)
    
    res
  })
  
  callModule(
    plotEIC_module, "peak_eic",
    react_x = reactive({
      req(peakSelected()$eic)
      peakSelected()$eic
    }),
    react_vline = reactive(
      peakSelected()$rtVline
    ),
    react_initXlim = reactive(
      peakSelected()$rtrange
    ) 
  )
}
