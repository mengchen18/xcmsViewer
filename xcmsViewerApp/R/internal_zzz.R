#' parse rt and mz range used by shiny app
#' @param x a character to be parsed, such as "123-133" or "133"
#' @param expand when x is parsed to a single value, how to expand the range
parseRange <- function(x, expand = 0.05) {
  v <- as.numeric(strsplit(x, "-")[[1]])
  if (any(is.na(v))) {
    warning('Unknown format x!')
    v <- c(0, 0)
  } else if (length(v) == 0) {
    v <- c(-Inf, Inf)
  } else if (length(v) == 2) {
    if (v[1] >= v[2]) {
      warning("Wrong range definition!")
      v <- c(0, 0)
    }
  } else if (length(v) == 1) {
    v <- c(max(0, v-expand), v+expand)
  } else {
    warning('Unknown format x!')
    v <- c(0, 0)
  }
  v
}

#' beeswarm or bar plot to show the intensity, used by shiny app
#' @param x numerical value
#' @param col color 
#' @param f value group
#' @param diffPheno pheno to distinguish

bs_or_bar <- function(x, col, f, diffPheno) {

  col <- col[diffPheno]  
  par(mar = c(4, 4, 1, 1))

  if (max(table(f)) > 1) {
    # beeswarm
    beeswarm::bxplot(x ~ f, col = "gray", ylab = "Intensity")
    beeswarm::beeswarm(x ~ f, pwcol = col, pch = 19, corral = "random", add = TRUE)
  } else {
    # barplot
    barplot(x, col = col, names.arg = f, las = 2, ylab = "Intensity")
  }
}
