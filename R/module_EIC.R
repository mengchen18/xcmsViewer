EIC_ui <- function(id, layoutFun = function(EICplot, colorControl, repPeakControl) {
  tagList(
    column(
      1,
      shinyWidgets::dropdownButton(
        circle = FALSE, status = "default", icon = icon("gear"), width = "800px", right = FALSE,
        tooltip = shinyWidgets::tooltipOptions(title = "Click to modify colors for EIC!"),
        br(),
        colorControl, # place holder 2 - color control widget,
        repPeakControl # place holder 3 - represetative control widget
      )
    ),
    column(
      11,
      EICplot # place holder 1 - EIC plot
    )
  )}
) {
  ns <- NS(id)
  EICplot <- plotlyOutput(ns("ms1eic"), height = "320px")
  colorControl <- ExpressionSetViewer::triselector_ui(ns("eic_color"))
  repPeakControl <- shinyWidgets::switchInput( inputId = ns("repeic"), size = "mini", label = "Representative EIC only", value = FALSE, labelWidth = 450)
  layoutFun(EICplot, colorControl, repPeakControl)
}

#' Module for EIC
#' @param input input
#' @param output output
#' @param session session
#' @param switchOnRepPeak Integer; if the numbe of sample is large, then only the representative
#'   EICs are shown. This parameter tells how many samples will be considered as "large size" of samples. 
#' @param pdata phenotype data
#' @param object xcmsViewer object
#' @param obj object returned by "xcmsAnnotationTab_module", the list object when a feature is selected,
#'   it should contain at list an element named "eic_param", which gives the range of mz (numeric vector of 
#'   length two; named mz) and range of RT (numeric vector of length two; named rt).
#' @param file integer; file index
EIC_module <- function(
  input, output, session, pdata, object, obj, switchOnRepPeak = 12, file = reactive(NULL)
) {
  
  observe({
    if (nrow(pdata()) > switchOnRepPeak)
      shinyWidgets::updateSwitchInput(session, inputId = "repeic", value = TRUE)
  }, priority = 1000)
  
  triset <- reactive({
    str_split_fixed(colnames(pdata()), "\\|", n = 3)
  })
  v1 <- callModule(
    ExpressionSetViewer::triselector_module, id = "eic_color", reactive_x = triset, label = "Color"
  )
  colvar <- reactiveVal("EIC")
  eiccol <- reactiveVal(c(EIC = "gray"))
  observe({
    cc <- paste(v1(), collapse = "|")
    if (!cc %in% colnames(pdata())) {
      colvar("EIC")
      eiccol(c(EIC = "gray"))
      return()
    }
    colvar(pdata()[, cc])
    cc <- sort(distinctColorPalette(k = length(unique(colvar()))))
    names(cc) <- unique(colvar())
    eiccol(cc)
  })
  
  pks <- reactive({
    req(obj()$featureid)
    pidx <- getFeatureMeta(object(), obj()$featureid)$"General|Extended|peakidx"[[1]]
    getPeak(object(), rowid = pidx)
  })
  
  eic <- reactive({
    req(pks <- pks())
    pks <- pks[pks$masterPeak == "+", ]
    file <- NA
    if (!is.null(file())) {
      file <- file()
    } else if (input$repeic) {
      ii <- c(
        which.max(pks$"General|All|rtmed"), which.min(pks$"General|All|rtmed"), which.max(pks$sn), 
        which.min(pks$sn), which.max(pks$into), which.min(pks$into))
      ii <- pks$sample[ii]
      for (ic in colvar()) {
        tu <- pks[colvar() == ic, ]
        ii <- c(ii, tu$sample[which.max(tu$into)])
      }
      file <- unique(na.omit(ii))
    }
    df <- getEICFromFeatureID(object(), ID = obj()$featureid)
    if (!is.na(file))
      df <- df[which(df$file %in% file), ]
    req(nrow(df) > 1)
    cv <- colvar()
    if (cv[1] != "EIC" || length(cv) > 1)
      cv <- cv[df$file]
    df$colorGroup <- cv
    df$tooltips <- pdata()[df$file, 1]
    df
  })
  
  output$ms1eic <- renderPlotly({
    req(eic())
    plotly_EIC(eic(), range = obj()$rtVline, legend.title = v1()[[3]], col = eiccol())
  })
  
  return(pks)
}