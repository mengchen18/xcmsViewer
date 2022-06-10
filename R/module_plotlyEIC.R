#' Plot EIC using plotly
#' @param x a data.frame with 5 columns
#'  - mz: mz
#'  - intensity: intensity
#'  - file: integer, the file index
#'  - colorGroup: character
#'  - tooltips: character tooltips
#' @param range a named numeric vector with three values named as "rtmin", "rtmax" and "rtmed"
#' @param xlim xlim
#' @param legend.title title of legend
#' @param col color code for EICs, it needs to be a named vector of colors, the names should be 
#'   the unique values of "x$colorGroup".
#' @importFrom randomcoloR distinctColorPalette
#' @rawNamespace import(plotly, except = groups)

plotly_EIC <- function(x, range, xlim = NULL, legend.title = "Group", col=NA) {
  
  im <- max(x$intensity)*1.05
  if (is.na(col[1])) {
    ucv <- unique(x$colorGroup)
    cc <- sort(distinctColorPalette(k = length(ucv)))
    names(cc) <- ucv
    col <- cc
  }
  col <- col[x$colorGroup]
  fig <- plot_ly()
  
  if (!is.null(range)) {
    rg <- c(range["rtmin"], range["rtmin"], range["rtmax"], range["rtmax"])
    if (all(!is.na(rg))) 
      fig <- add_trace(
        fig,
        x = rg,
        y = c(0, im, im, 0),
        type = 'scatter',
        mode = "lines",
        fill = 'toself',
        fillcolor = "RT range",
        hoveron = 'points+fills',
        line = list(
          color = rgb(202, 255, 112, 75, maxColorValue = 255)
        ),
        showlegend = FALSE,
        text = "RT range",
        hoverinfo = 'text'
      )
    
    if (!is.na(range["rtmed"])) 
      fig <- add_trace(
        fig,
        x = c(range["rtmed"], range["rtmed"]),
        y = c(0, im),
        line = list(
          color = rgb(150, 150, 150, maxColorValue = 255)
        ),
        type = "scatter",
        mode = "lines",
        name = "median RT",
        showlegend = FALSE)
  }
  
  i0 <- rep(TRUE, nrow(x))
  
  if (length(unique(x$colorGroup)) > 1) {
    for (ii in unique(x$colorGroup)) {
      i3 <- which(x$colorGroup == ii)
      i4 <- which(x$file[i3] == x$file[i3[1]])
      it <- i3[i4]
      i0[it] <- FALSE
      fig <- add_trace(fig, data = x[it, ], x = ~ rt, y = ~ intensity, 
                       line = list(color = col[it[1]]),
                       type = 'scatter', mode = 'lines', name = ii,
                       showlegend = TRUE, hoveron = 'points+lines',
                       text = ~ paste(colorGroup, tooltips, sep = "\n")
      )
    }
  }
  
  if(any(i0))
    fig <- add_trace(fig, data = x[i0, ], x = ~ rt, y = ~ intensity,
                     color = ~ I(col[i0]), split = ~ file,
                     type = 'scatter', mode = 'lines',
                     showlegend = FALSE, hoveron = 'points+lines',
                     text = ~ paste(colorGroup, tooltips, sep = "\n")
    )
  
  fig <- layout(
    fig, 
    legend=list(
      title=list(text=sprintf('<b>%s </b>', legend.title)),
      orientation="h",
      yanchor="bottom",
      y=1.02,
      xanchor="right",
      x=1
    ),
    xaxis = list(title = "Retention time", range = xlim)
  )
  fig
}
