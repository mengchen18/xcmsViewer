xcmsAnnotationTab_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # select features
    absolutePanel(
      top = 392, right = 15, style = "z-index: 999;",
      div(style="display: inline-block;vertical-align:top;",
        downloadButton(ns("downloadMS2"), label = "Spectra (MGF)")),
      div(style="display: inline-block;vertical-align:top;",
        downloadButton(ns("downloadAllMS2"), label = "All MS2 Spectra (MGF)")
        )
      ),    
    fluidRow(
      EIC_ui(ns("eic0")),
      HTML("<br>"),
      tabsetPanel(
        tabPanel(
          "MS2 & annotation",
          msAnnotation_UI(ns("ms2v"))
        ),
        tabPanel(
          "Peaks",
          DT::dataTableOutput(ns("peaktable")),
          plotlyOutput(ns("ms1eic_indpeak"), height = "320px")
        )
      )
    )    
  )
}

xcmsAnnotationTab_module <- function(
  input, output, session, pdata, fdata, expr, feature_selected, sample_selected, object
) {
  
  ns <- session$ns  

  .object <- reactive({
    req(object())
    object()
    })
  
  ci <- c("ID", "rtmed", "mzmed", "annot_ms1", "annot_ms2", "npeaks", "purity", "rtmin", "rtmax", "mzmin", "mzmax")
  ftab <- reactive({
    req(length(feature_selected()) == 1 && !is.logical(feature_selected()))
    tab <- fdata()[feature_selected(), ]
    colnames(tab) <- sapply(strsplit(colnames(tab), "\\|"), "[", 3)
    tab <- tab[, c(ci, "ms2spectrum")]
    ic <- sapply(tab, is.numeric)
    tab[ic] <- lapply(tab[ic], round, digits = 3)
    tab
  })
  
  # return items 
  obj <- eventReactive(ftab(), {

    req(length(feature_selected()) == 1 && !is.logical(feature_selected()))

    f <- ftab()
    res <- list()
    
    # 1. row names of feature table (feature name)
    res$featureid <- f$ID # feature_selected()
    if (is.na(res$featureid))
      req(NULL)
    # print(sprintf("feature ID selected: %s", res$featureid))
    
    # 3. vline in eic plot
    res$rtVline <- c(
      rtmed = f$rtmed,
      rtmin = f$rtmin,
      rtmax = f$rtmax)
    
    # 4. peak table
    v <- getPeak(.object(), rowid = f$peakidx[[1]])
    v[sapply(v, is.numeric)] <- lapply(v[sapply(v, is.numeric)], round, digits = 4)
    res$peakTab <- v
    # 5. annotation
    res$annot <- getAnnot(.object(), ID = res$featureid)
    # 6. MS2 scan and meta
    sids <- getScanIDFromFeatureID(object = .object(), ID = res$featureid)
    res$ms2scanMeta <- getScanMeta(.object(), sids)
    res$ms2scan <- getScan2(object = .object(), scanId = sids)
    # 7. consensus spectrum
    res$consensusMS2Spectrum <- str2spectra(f$ms2spectrum)
    res
  })
  
  # =============== download MS2 MGF ==============
  mod <- reactive({
    req(ad <- getAnnot(.object())$Adduct)
    nad <- nchar(ad[1])
    paste0(1, substr(ad[1], nad, nad))
  })
  mgf <- reactive({
    prepareMGF(scan = obj()$ms2scan, scanMeta = obj()$ms2scanMeta, 
               cons = obj()$consensusMS2Spectrum, mode = mod(), 
               featureId = obj()$featureid)
  })
  output$downloadMS2 <- downloadHandler(
    filename = function() {
      paste('MGF', Sys.Date(), '.mgf', sep='')
    },
    content = function(con) {
      if (length(mgf()) == 0)
        return(NULL)
      writeLines(mgf(), con)
    },
    contentType = "text/mgf"
  )

   dataModal <- function(empty = FALSE) {
      modalDialog(
        title = ifelse(empty, 
          "Nothing to export!",
          "Exporting MGFs ..."),        
        if (!empty)
          "This may take a a few minutes ...",
        easyClose = empty,
        footer = NULL
      )
    }

  output$downloadAllMS2 <- downloadHandler(
    filename = function() {
      paste('allMS2Spec', Sys.Date(), '.mgf', sep='')
    },
    content = function(con) {
      showModal(dataModal())
      q <- prepareAllMGFs(.object(), con)
      if (file.exists(con)) 
        removeModal() else
          showModal(dataModal(empty = TRUE))
    },
    contentType = "text/mgf"
  )
  
  # ================ MS1 EIC plot ===================
  pks <- callModule(
    EIC_module, id = "eic0", pdata = pdata, object = .object, obj = obj, switchOnRepPeak = 12
    )
  ## =============== MS1 peaks ==================
  output$peaktable <- DT::renderDataTable(
    formatFeaturePeakTab(pks())
    )

  eic_indpeak <- reactive({
    req(i <- input$peaktable_rows_selected)
    eic <- getEICFromFeatureID(.object(), ID = obj()$featureid)
    p1 <- pks()[i, ]
    dd <- eic[eic$file == p1$sample, ]
    dd$colorGroup <- "EIC"
    dd$tooltips <- pdata()[dd$file, 1]
    rtr <- c(min(pks()$rtmin)-10, max(pks()$rtmax)+10)
    print(rtr)
    list(peaks = dd, 
         vline = c("rtmin" = p1$rtmin, 
                   "rtmax" = p1$rtmax,
                   "rtmed" = p1$rt),
         xlim = rtr)
  })
  
  output$ms1eic_indpeak <- renderPlotly({
    req(eic_indpeak())
    plotly_EIC(eic_indpeak()$peaks, range = eic_indpeak()$vline, xlim = eic_indpeak()$xlim)
  })
  
  ## =============== MS2 visual =================
  callModule(msAnnotation, "ms2v", dat = .object, featureSelected = obj)
}
