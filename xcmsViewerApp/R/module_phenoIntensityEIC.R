phenoIntensityEIC_UI <- function(id) {
  ns <- NS(id)

  shiny::tagList(
    wellPanel(
          tags$div(
            style="margin-bottom:8px;",
            DT::dataTableOutput(ns("pheno")) ###
          ),
          fluidRow(
            column(
              6,
              checkboxInput(ns("multiSelect"), "Mutiple selection", value = FALSE)
            ),
            column(
              6, align = 'right',
              actionButton(ns("clearRow"), "Clear selection")
            )
          ),
          style = "font-size: 12px"
        ),
    wellPanel(
          selectizeInput(ns("phenoCols"), "Phenotype:", choices = NULL),
          plotOutput(ns("intensityPlot"), height = 300, click = ns("bs_or_bar_click")),
          plotEIC_UI(ns("eic_plot"))
        )
    )
}


phenoIntensityEIC <- function(input, output, session, dat, featureSelected = reactive(NULL)) {

  # 1. pheontype table
  output$pheno <- DT::renderDataTable({
    if (input$multiSelect)
    ss <- "multiple" else
      ss <- "single"
    DT::datatable(
      dat()$pheno,
      selection = ss,
      filter = "top",
      rownames = FALSE,
      class="table-bordered compact",
      caption = "Phenotype/Meta data",
      option = list(scrollX = TRUE, scrollY = "220px", dom = "t", paging = FALSE)
    )
  })
  dtProxy <- DT::dataTableProxy('pheno')
  observeEvent(input$clearRow, {
    DT::selectRows(dtProxy, NULL)
  })
  
  # 2. update phenotype of interest
  observe(
    updateSelectizeInput(session, "phenoCols", choices = colnames(dat()$pheno))
    )

  # 3. EIC plot
  qt <- reactive({
    req(input$phenoCols)
    as.character(dat()$pheno[[input$phenoCols]])
    })

  colvar <- callModule(
    plotEIC_module, "eic_plot",
    react_x = reactive( featureSelected()$eic ),
    react_diffPheno = qt,
    react_select = reactive( input$pheno_rows_selected ),
    react_vline = reactive( featureSelected()$rtVline ),
    reset_click = reactive( c(input$bs_or_bar_click, input$pheno_rows_selected) )
  )
  
  observeEvent(colvar()$file_clicked, {
    req(colvar()$file_clicked)
    DT::selectRows(dtProxy, colvar()$file_clicked)
  })  

  # 4. barplot/beeswarm
  bs_or_bar_data <- reactive({

    req(featureSelected()$featureid)

    ir <- fmatch(featureSelected()$featureid, dat()$features$intensities$ID)
    it <- setdiff(colnames(dat()$features$intensities), "ID")
    int <- unlist(dat()$features$intensities[ir, it])
    int[is.na(int)] <- 0

    ## there might be a bug
    fc <- dat()$pheno[[input$phenoCols]]
    if (length(fc) == 0 ) {
      fc <- rep("N/A", length(int))
    } else if (length(fc) != length(int)) {
      warning("check module_phenoIntensityEIC phenoCols problem!")
      fc <- rep("N/A", length(int))
    }

    df <- data.frame(
      i = int, 
      f = fc,
      stringsAsFactors = FALSE
      )

    if (max(table(df$f)) > 1) {
      # beeswarm
      coord <- beeswarm(df$i ~ df$f, do.plot = FALSE)
      rvec <- paste(coord$x.orig, coord$y.orig)
      coord$ID <- fmatch(rvec, paste(df$f, df$i))
    } else {
      xc <- barplot(df$i)
      coord <- data.frame(
        x = rep(xc, each = 20),
        y = unlist(lapply(df$i, seq, 0, length.out = 20)),
        ID = rep(1:length(xc), each = 20)
        )
    }
    list(df=df, coord=coord)
    })

  observe({
    cat("Phenotype to distinguish:\n")
    print(colvar()$diffPheno)
    })

  output$intensityPlot <- renderPlot({
    bs_or_bar(
      bs_or_bar_data()$df$i, 
      col = colvar()$col,
      f = bs_or_bar_data()$df$f,
      diffPheno = colvar()$diffPheno
    )
  })
  
  observe({
    req(input$bs_or_bar_click)
    np <- nearPoints(
      bs_or_bar_data()$coord, input$bs_or_bar_click, xvar = "x", yvar = "y", maxpoints = 1, threshold = 25
      )
    DT::selectRows(dtProxy, np$ID)
  })
}