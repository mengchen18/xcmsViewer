#' Check if an XIC looks like a Guassian shape
#' @param x a data.frame with two columns named as 'mz' and 'int' (for intensity)
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom stats deviance dist median na.omit sd var
#' @importFrom Biobase pData
#'
.evalPeak <- function(x) {
  
  nlsrsq <- function (mdl, y, param) {
    adj <- (sum(!is.na(y)) - 1)/(sum(!is.na(y)) - param)
    sum.sq <- (sum(!is.na(y)) - 1) * var(y, na.rm = TRUE)
    rsq <- 1 - (adj * (deviance(mdl)/sum.sq))
    return(rsq)
  }
  
  tab <- na.omit(x)
  colnames(tab)[colnames(tab) == "intensity"] <- "int"
  
  if (nrow(tab) < 5) {
    res <- NA 
  } else {
    res <- try(
      minpack.lm::nlsLM(
        int ~ k*exp(-((rt-mu)^2)/(2*sigma^2)) / sqrt(2*pi*sigma^2)+b, 
        start=c(mu=mean(tab$rt),sigma=sd(tab$rt), k = max(tab$int), b = 0), 
        data = tab, control = nls.lm.control(maxiter = 500)
      ), silent = TRUE)
    if (inherits(res, "try-error")) 
      res <- NA
  }
  
  rsq <- rtgap <- intgap <- rtintgap  <- truncated <- b <- NA
  
  if (!is.na(res)[1]) {
    rsq <- nlsrsq(res, tab$int, param  = 4)
    
    d <- dist(tab$rt)
    rtgap <- max(d)/median(d)
    
    d <- dist(tab$int)
    intgap <- max(d)/median(d)
    
    tabn <- apply(tab, 2, function(x) {
      x <- x-min(x)
      x/max(x)
    })
    d <- dist(tabn)
    rtintgap <- max(d)/median(d)
    
    estmu <- summary(res)$parameters["mu", "Estimate"]
    if (estmu < min(tab$rt) || estmu > max(tab$rt)) 
      truncated <- Inf else
        truncated <- log10((max(tab$rt) - estmu)/(estmu - min(tab$rt)))
    
    b <- summary(res)$parameters["b", "Estimate"]
  }
  
  list(
    model = res,
    stats = c(
      "rsq" = rsq,
      "rtgap" = rtgap,
      "intgap" = intgap,
      "rtintgap" = rtintgap,
      "truncated" = truncated,
      "b" = b
    )
  )  
}
