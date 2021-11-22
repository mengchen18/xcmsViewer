msAnnotation_UI <- function(id) {
  
  ns <- NS(id)
  shinyjs::useShinyjs()
  
  shiny::tagList(
    fluidRow(
      column(
        6,
        wellPanel(
          DT::dataTableOutput(ns("massAnnotTab")),
          mirrorPlotUI(ns("ms2mirrorplot")),
          style = "padding: 5px 15px 5px 15px; font-size: 12px"
        )
      ),
      column(
        6,         
        wellPanel(
          DT::dataTableOutput(ns("ms2table_scans")),
          mirrorPlotUI(ns("ms2scanPeaksPlot")),
          style = "padding: 5px 15px 5px 15px; font-size: 12px"
        )
      )
    )
  )
}

msAnnotation <- function(input, output, session, dat, featureSelected=reactive(NULL)) {
  
  ## ================= 1. MS1 annotation table =================
  maTab <- reactive({
    an <- featureSelected()$annot
    an$ms1
  })
  
  createLink <- function(val) {
    v <- sprintf(
      '<a href="https://pubchem.ncbi.nlm.nih.gov/compound/%s" target="_blank" class="btn btn-primary">%s</a>',
      val, val)
    v[is.na(val)] <- NA
    v
  }
  
  output$massAnnotTab <- DT::renderDataTable({
    req(maTab())
    tab <- maTab()[, setdiff(colnames(maTab()), "ID")]
    tab$CID <- createLink(tab$CID)
    dt <- DT::datatable(    
      tab,
      selection = list(mode="single", selected = 1),
      rownames = FALSE,
      filter = "top",
      class="table-bordered compact",
      options = list(scrollX = TRUE, dom = 't', scrollY = "180px", paging = FALSE),
      escape=FALSE,
      caption = sprintf("Annotation using exact mass matching, feature %s", featureSelected()$featureid)
    )
    dt
  }, escape=FALSE)
  
  ## ================= 2. MS2 scan table ================= 
  # scan meta table
  ms2scans <- reactive({
    t <- featureSelected()$ms2scanMeta[, c("ID", "rt", "tic", "peakCount", "fromFile", "precMz")]
    t[t$ID %in% featureSelected()$ms2scan$ID, ]
  })
  
  formatMs2Tab <- function(tab) {
    dt <- DT::datatable( 
      tab,
      selection =  list(mode="single", selected = 1),
      rownames = FALSE,
      filter = "top",
      class="table-bordered compact",
      caption = sprintf("MS2 of feature %s", featureSelected()$featureid),
      options = list(scrollX = TRUE, scrollY = "180px", paging = FALSE, dom = 't')
    )
    DT::formatStyle(dt, columns = 1:ncol(tab), fontSize = '85%')
  }
  
  output$ms2table_scans <- DT::renderDataTable({
    formatMs2Tab(ms2scans())
  })
  
  ms2scanPeaks <- reactive({
    req(i <- input$ms2table_scans_rows_selected)
    req(scn <- featureSelected()$ms2scan)
    f <- ms2scans()[i, ]
    m_peaks <- scn[scn$ID == f$ID, ]
    m_text <- sprintf(
      "Prec MZ: %s,RT:%s; TIC:%s; File:%s, Peak count: %s",
      round(f$precMz, 4), 
      round(f$rt, 4), 
      f$tic, f$fromFile, f$peakCount)
    list(measured =  m_peaks,
         measuredLegend  = m_text)
  })
  
  callModule(mirrorPlotModule, "ms2scanPeaksPlot",
             measured = reactive(ms2scanPeaks()$measured),
             standard = reactive(NULL),
             legend.measured  = reactive(ms2scanPeaks()$measuredLegend),
             legend.standard = reactive(NULL),
             ppmtol = reactive(getMassTol(dat())["MS2"]))
  
  ## ================= 3. MS2 scan table - linked to mirror plot =================
  # reference peak list
  mirrorData <- reactive({
    req(i <- input$massAnnotTab_rows_selected)
    t <- maTab()[i, ]
    ref <- featureSelected()$annot$ms2[[t$InChIKey]]
    ref <- na.omit(ref)
    # show consensus spectrum if no 
    if (is.null(ref) || nrow(ref) == 0)
      ref <- NULL
    list(
      ref = ref,
      measured = featureSelected()$consensusMS2Spectrum
    )
  })
  
  callModule(mirrorPlotModule, "ms2mirrorplot",
             measured = reactive(mirrorData()$measured),
             standard = reactive(mirrorData()$ref),
             legend.measured  = reactive("Consensus MS2 spectrum of measured"),
             legend.standard = reactive("Reference consensus spectrum"),
             ppmtol = reactive(getMassTol(dat())["MS2"]))
}

