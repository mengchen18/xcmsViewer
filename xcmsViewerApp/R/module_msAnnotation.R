msAnnotation_UI <- function(id) {
	ns <- NS(id)
	shiny::tagList(
		wellPanel(
              DT::dataTableOutput(ns("massAnnotTab")),
              style = "padding: 5px 15px 5px 15px; font-size: 12px"
            ),
            wellPanel(
              DT::dataTableOutput(ns("ms2table")),
              style = "padding: 5px 15px 5px 15px; font-size: 12px"
            ),
            wellPanel(
              mirrorPlotUI(ns("ms2mirrorplot")))
		)
}

msAnnotation <- function(input, output, session, dat, featureSelected=reactive(NULL)) {

  maTab <- reactive({
    i_tab <- data.frame(
      Name = character(),
      metaID = character(),
      Score = numeric(),
      MS2cos = numeric(),
      Formula = character(),
      Adduct = character(),
      Delta = numeric(),
      stringsAsFactors = FALSE
    )
    
    if (!is.null(featureSelected()$ms1annot)) {
      if (nrow(featureSelected()$ms1annot) > 0) {
        i_tab <- rbind(
          i_tab, 
          data.frame(
            Name = featureSelected()$ms1annot$name,
            metaID = featureSelected()$ms1annot$accession,
            Score = round(featureSelected()$ms1annot$atScore, digits=3),
            MS2cos = round(featureSelected()$ms1annot$cos, digits = 3),
            Formula = featureSelected()$ms1annot$chemical_formula,
            Adduct = featureSelected()$ms1annot$adduct,
            Delta = round(featureSelected()$ms1annot$deltaPPM, digits = 3),
            stringsAsFactors = FALSE
          )
        )
      }
    }
    i_tab
    })

	output$massAnnotTab <- DT::renderDataTable({
    v <- NULL
    if (nrow(maTab()) > 0)
      v <- 1

    dt <- DT::datatable(    
      maTab(),
      selection = list(mode="single", selected = v),
      rownames = FALSE,
      filter = "top",
      class="table-bordered compact",
      options = list(scrollX = TRUE, dom = 't', scrollY = "200px", paging = FALSE),
      caption = "Annotation using exact mass matching"
    )
    dt
  })
  maTabProxy <- DT::dataTableProxy('massAnnotTab')

  observeEvent(maTab(), {
    req(maTab())
    req(ms2tab())
    req(i <- input$massAnnotTab_rows_selected)
    req(id <- maTab()$metaID[i])
    ii <- fastmatch::fmatch(id, ms2tab()$metaID)
    if (is.na(ii)) req(NULL)
    DT::selectRows(ms2TabProxy, ii)
    })

  observeEvent(input$massAnnotTab_cell_clicked, {    
    req(maTab())
    req(ms2tab())
    req(i <- input$massAnnotTab_rows_selected)
    req(id <- maTab()$metaID[i])
    ii <- fastmatch::fmatch(id, ms2tab()$metaID)
    if (is.na(ii)) {
      ii <- NULL
    } else if (!is.null(input$ms2table_rows_selected)) {
      idms2 <- ms2tab()$metaID[input$ms2table_rows_selected]
      if (is.na(idms2)) {
        ii <- NULL
      } else if (id == idms2) {
        ii <- input$ms2table_rows_selected
      }      
    }
    DT::selectRows(ms2TabProxy, ii)
    })

  
  # 2. MS2 scan table
  ms2tab <- reactive({

    ms2tab <- featureSelected()$ms2scan
    
    i_tab <- data.frame(
      metaID = character(),
      ms2ID = character(),
      cos = numeric(),
      refID = character(),
      stringsAsFactors = FALSE
    )
    if (!is.null(ms2tab)) {
      if (nrow( ms2tab > 0)) {
        i_tab <- rbind(
          i_tab,
          data.frame(
            metaID = ms2tab$database_id,
            ms2ID = ms2tab$query_peaks,
            cos = round(ms2tab$cos, digits = 3),
            refID = ms2tab$annot_peaks,
            stringsAsFactors = FALSE
          )
        )
      }
    }
    i_tab
  })

  output$ms2table <- DT::renderDataTable( {    

    req(ms2tab())
    v <- NULL
    if (nrow(ms2tab()) > 0)
      v <- 1
    dt <- DT::datatable( 
      ms2tab(),
      selection =  list(mode="single", selected = v),
      rownames = FALSE,
      filter = "top",
      class="table-bordered compact",
      options = list(scrollX = TRUE, dom = 't', scrollY = "180px", paging = FALSE),
      caption = "MS2 related to this feature"
    )
    DT::formatStyle(dt, columns = 1:ncol(ms2tab()), fontSize = '85%')
  })

  ms2TabProxy <- DT::dataTableProxy('ms2table')


  # observeEvent(input$ms2table_rows_selected, {
  observeEvent(input$ms2table_cell_clicked, {
    req(maTab())
    req(ms2tab())
    req(i <- input$ms2table_rows_selected)
    id <- ms2tab()$metaID[i]
    if (is.na(id))
      ss <- NULL else
        ss <- fastmatch::fmatch(id, maTab()$metaID)
    DT::selectRows(maTabProxy, ss)
  })
  
  # 3. MS2/mirror plot
  ms2peaks <- reactive({
    req( i <- input$ms2table_rows_selected )
    f <- featureSelected()$ms2scan[i, ]
    req(f)
    if (is.na(f$query_peaks))
      req(NULL)
    
    i1 <- dat()$scanIntensityTab$ID == f$query_peaks
    m_peaks <- dat()$scanIntensityTab[i1, ]
    if (!is.na(f$database_id)) {
      s_peaks <- dat()$matchedRefFragments$peakList[[f$annot_peaks]]
      i3 <- fmatch(f$annot_peaks, dat()$matchedRefFragments$meta$id2)
      s_meta <- dat()$matchedRefFragments$meta[i3, ]
      s_text <- sprintf("%s; CE-level: %s; i-mode: %s; ID: %s",
                        s_meta$notes, s_meta$"collision_energy_level", s_meta$"ionization_mode",
                        s_meta$"database_id")
    } else {
      s_peaks <- NULL
      s_text <- ''
    }
    
    iq <- grep(f$query_peaks, featureSelected()$peakTab$ms2Scan)
    precmz <- featureSelected()$peakTab$mz[iq]
    scanmet <- dat()$scanMetaTab[fmatch(f$query_peaks, dat()$scanMetaTab$ID), ]
    m_text <- sprintf("Prec MZ: %s,RT:%s; TIC:%s; File:%s, Peak count: %s", 
                      precmz, scanmet$rt, scanmet$tic, scanmet$fromFile, scanmet$peakCount)
    
    list(measured = m_peaks, standard = s_peaks, standardLegend = s_text, measuredLegend = m_text)
  })
  
  callModule(mirrorPlotModule, "ms2mirrorplot",
             measured = reactive(ms2peaks()$measured),
             standard = reactive(ms2peaks()$standard),
             legend.measured  = reactive(ms2peaks()$measuredLegend),
             legend.standard = reactive(ms2peaks()$standardLegend),
             ppmtol = reactive(10))
}
