featureStatsTab_UI <- function(id) {
  ns<- NS(id)
  wellPanel(
    tabsetPanel(
      tabPanel(
        "Filter table",
        fluidRow(
            column(
              3, textInput(ns("range_RT"), label = "Range RT", placeholder = "e.g. 120-130")
            ),
            column(
              3, textInput(ns("range_MZ"), label = "Range MZ", placeholder = "e.g. 320-321 or 455")
            ),
            column(
              4, textInput(ns("search_annot"), label = "Search metabolite", placeholder = "e.g. proline or prol")
            ),
            column(
              2, tags$div(
                style="margin-top:25px;",
                actionButton(ns("sbutton"), "Filter")
                )
              )
          )
        ),
      tabPanel(
        "Customize table",
        selectInput(ns("displayCols"), "Select columns to show", choices = NULL, multiple = TRUE)
        )
      ),
    DT::dataTableOutput(ns("featureTable"))
    )
}


featureStatsTab <- function(input, output, session, dat, dataChanged) {

  featureTab <- reactive({
    req(dat())

    input$sbutton
      isolate({
        r_rt <- parseRange(input$range_RT)
        r_mz <- parseRange(input$range_MZ)
        s_meta <- input$search_annot    
        i1 <- dat()$features$meta$mzmed > r_mz[1] & dat()$features$meta$mzmed < r_mz[2]
        i2 <- dat()$features$meta$rtmed > r_rt[1] & dat()$features$meta$rtmed < r_rt[2]
        i3 <- TRUE
        if (s_meta != "")
          i3 <- sapply(dat()$annotationMass, function(x) {
            any(grepl(s_meta, x$name, ignore.case = TRUE))
          })
      })
    dat()$features$meta[i1 & i2 & i3, ]
  })

  observe({
    updateSelectInput(session, "displayCols", choices = colnames(featureTab()), 
      selected = intersect(
        c("ID", "Annotation", "QC", "mzmed", "rtmed"), 
        colnames(featureTab())
        )
      )
  })
  
  ## render the table
  output$featureTable <- DT::renderDataTable({
    req( input$displayCols )
    ii <- input$displayCols %in% colnames(featureTab())
    req( ic <-  input$displayCols[ii] )
    dt <- DT::datatable(
      featureTab()[, ic],
      selection = "single",
      rownames = FALSE,
      filter = "top",
      options = list(scrollX = TRUE, dom = 'tip', pageLength=25)
    )
    dt <- DT::formatRound(dt, columns = 4:7, digits = 4)
    DT::formatStyle(dt, columns = 2, fontSize = '85%')
  })
  
  featureRowSelected <- reactiveVal(NULL)
  observeEvent(list(input$sbutton, dataChanged(), featureTab(), dat()), {
    print("feature selected set to NULL!")
    featureRowSelected(NULL) 
    })
  observe( 
    featureRowSelected( input$featureTable_rows_selected )
    )
  
  ## after select a row from the table, returns
  reactive({

    req( featureRowSelected() )
    res <- list()
    
    f <- featureTab()[featureRowSelected(), ]
    # 1. row names of feature table (feature name)
    res$featureid <- f$ID
    if (is.na(res$featureid))
      req(NULL)
    print(sprintf("feature ID selected: %s", res$featureid))
    
    # 2. feature EIC
    mzoff <- 25*1e-6*f$mzmed
    res$eic <- eic(
      itab = dat()$scanIntensityTab, 
      mtab = dat()$scanMetaTab, 
      rt = c(f$rtmin - 60, f$rtmax + 60), 
      mz = c(f$mzmin - mzoff, f$mzmax + mzoff)
      )
    
    # 3. vline in eic plot
    res$rtVline <- c(
      rtmed = f$rtmed,
      rtmin = f$rtmin,
      rtmax = f$rtmax)
    
    # 4. peak table
    v <- dat()$peaks[f$peakidx[[1]], ]
    v[sapply(v, is.numeric)] <- lapply(v[sapply(v, is.numeric)], round, digits = 4)
    res$peakTab <- v
    
    # 5. MS1 annotation table
    res$ms1annot <- dat()$annotationMass[[res$featureid]]
    
    # 6. MS2 scans 
    ms2scan <- unlist(strsplit(res$peakTab$ms2Scan, ";"))
    ms2scan <- ms2scan[ms2scan %fin% dat()$scanIntensityTab$ID]
    if (length(ms2scan) == 0)
      matched <- NULL else {
        i2 <- dat()$annotationFragment$query_peaks %fin% ms2scan
        matched <- dat()$annotationFragment[i2, c("query_peaks", "annot_peaks", "cos", "database_id")]
        v <- setdiff(ms2scan, matched$query_peaks)
        if (length(v) > 0) {
          ms2scandf <- data.frame(query_peaks = v, 
                                  annot_peaks = NA, 
                                  cos = NA, 
                                  database_id = NA,
                                  stringsAsFactors = FALSE, 
                                  check.names = FALSE)
          matched <- rbind(matched, ms2scandf)
        }
      } 
    res$ms2scan <- matched
    res
  })
}