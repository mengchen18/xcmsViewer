#' Functions for multiple t-test
#' @param x the expression matrix, should be log10 transformed
#' @param pheno the phenotype data.frame
#' @param compare comparisons to be done. An nx3 matrix, the first colume indicate the 
#'   phenotype column, the 2nd and 3rd columns indicate the two groups that should be compared
#' @param median.center logical; whether the expression matrix should be median centered
#' @param fillNA whether the NA should be filled
#' @param ... other parameters passed to t.test
#' @importFrom omicsViewer fillNA
#' @importFrom matrixStats colMedians
#' 
multi.t.test2 <- function(x, pheno, compare = NULL, median.center = FALSE, fillNA = TRUE, ...) {
  
  if (ncol(x) == 1)
    fillNA <- FALSE
  x.raw <- x
  if (is.null(rownames(x))) rownames(x) <- 1:nrow(x)
  if (fillNA) x <- fillNA(x)
  if (median.center)
    x <- sweep(x, 2, colMedians(x, na.rm = TRUE), "-") + median(x, na.rm = TRUE)
  
  if (is.null(compare)) {
    return( list(ttest = NULL, mat = x) )
  } else {
    tl <- lapply(unique(compare[, 1]), function(x) {
      x <- compare[x == compare[, 1], -1, drop = FALSE]
      unique(unlist(split(x, row(x))))
    })
    names(tl) <- unique(compare[, 1])
  }
  
  df <- data.frame(metabolite = rownames(x), stringsAsFactors = FALSE)
  for ( i in names(tl) ) {
    for ( j in tl[[i]] ) {
      m <- x[, pheno[[i]] == j, drop = FALSE]
      rv <- rowSums(!is.na(x.raw[, pheno[[i]] == j, drop = FALSE]))
      rm <- rowMeans(m, na.rm = TRUE)
      rq <- rank(rm, na.last = TRUE)/sum(!is.na(rm))
      rq[is.na(rm)] <- NA
      df[[paste("Stats|Mean", j, sep = "|")]] <- rm
      df[[paste("Stats|N value", j, sep = "|")]] <- rv
      df[[paste("Stats|Quantile", j, sep = "|")]] <- rq
    }
  }  
  for ( i in 1:nrow(compare) ) {
    v <- compare[i, ]
    i1 <- pheno[[v[1]]] == v[2]
    i2 <- pheno[[v[1]]] == v[3]
    
    tv <- apply(x, 1, function(xx) {
      # ve <- sum(!is.na(xx[i1])) < 2 || sum(!is.na(xx[i2])) < 2
      ve <- TRUE
      t <- try(t.test(xx[i1], xx[i2], var.equal = ve, ...), silent = TRUE)
      if (class(t) != "htest")
        return(c(pvalue = NA, mean.diff = NA))
      if (length(t$estimate) == 1)
        md <- t$estimate[[1]] else
          md <- t$estimate[1] - t$estimate[2]
        c(pvalue = t$p.value, mean.diff = md)
    })
    
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "pvalue", sep = "|")]] <- tv[1, ]
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "log.pvalue", sep = "|")]] <- -log10(tv[1, ])
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "fdr", sep = "|")]] <- p.adjust(tv[1, ], method = "fdr")
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "log.fdr", sep = "|")]] <- -log10(df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "fdr", sep = "|")]])
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "mean.diff", sep = "|")]] <- tv[2, ]
  }  
  list(ttest = df, mat = x)
}

#####
#' Function to perform PCA
#' @description the PCA result is ready to be included in the prunedXcmsSet
#' @param x the expression matrix
#' @param n number of components to include
#' @param prefix the prefix of headers in loading and PC matrix
#' @param fillNA logical
#' @param samples which samples to be included in the PCA
#' 
pca <- function(x, n = 6, prefix = "", fillNA = TRUE, samples = NULL) {
  
  if (!is.null(samples)) {
    if (is.logical(samples))
      n <- min(length(samples), ncol(x), nrow(x)) else
        n <- length(samples)
  }
  
  temp_sample <- matrix(NA, nrow = ncol(x), ncol = n)  
  temp_feature <- matrix(NA, nrow = nrow(x), ncol = n)

  if (ncol(x) == 1)
    fillNA <- FALSE

  if (fillNA) x <- fillNA(x)
  if (!is.null(samples))
    x <- x[, samples, drop = FALSE] else
      samples <- 1:ncol(x)
  if (ncol(x) <= 2) 
    return(NULL)

  ir <- which(rowVars(x, na.rm = TRUE) > 0 & rowSums(is.na(x)) == 0)
  if (length(ir) <= 2)
    return(NULL)  

  writePC <- function(x, n) {
    var <- round(x$sdev[1:n]^2/(sum(x$sdev^2))*100, digits = 1)
    xx <- x$x[, 1:min(n, ncol(x$x))]
    colnames(xx) <- paste0(prefix, "|", colnames(xx), " (", var, "%", ")")
    pp <- x$rotation[, 1:min(n, ncol(x$x))]
    colnames(pp) <- paste0(prefix, "|", colnames(pp), " (", var, "%", ")")
    list(samples = xx, features = pp)
  }
  
  pc <- prcomp(t(x[ir, ]))
  pc <- writePC(pc, n = n)
  
  temp_sample[ samples, ] <- pc$samples  
  colnames(temp_sample) <- colnames(pc$samples)
  temp_feature[ir, ] <- pc$features
  colnames(temp_feature) <- colnames(pc$features)
  list(samples = temp_sample, features = temp_feature) 
}

#' Function to perform multiple PCA
#' @description the PCA result is ready to be included in the prunedXcmsSet
#' @param x the expression matrix
#' @param pheno the phenotype data.frame
#' @param compare A named list object tells how the PCA should be performed. 
#'   for example:
#'   list(subsetName1 = c("phenotypeHeader1", "value1", "value2", "value3"),
#'        subsetName2 = c("phenotypeHeader2", "cat1", "cat2"))
#'   indicate the samples in "phenotypeHeader1" column of phenotype mapped to 
#'   value1, value2 and value3 will be included in the PCA. Similarly, the 
#'   samples in "phenotypeHeader2" column of phenotype mapped to 
#'   cat1 and cat2 will be included in the PCA, as a separate PCA.
#' @param n number of components to include
#' @param fillNA logical
#' @param prefix the prefix used in the header for columns of PCA results
#' @importFrom matrixStats rowVars
#' 
multi.pca <- function(x, pheno, compare = list("_all_" = TRUE), n = 6, fillNA = TRUE, prefix = "PCA") {
  
  cn <- setdiff(setdiff(sapply(compare, "[", 1), TRUE), c("_all_", colnames(pheno)))
  if (length(cn) > 0)
    stop(sprintf("These columns do not exist in pheno: %s", paste(cn, collapse = ", ")))
  
  if (!inherits(x, "matrix"))
    x <- apply(x, 2, as.numeric)
  
  t <- lapply(names(compare), function(i) {
    if (i == "_all_") {
      v <- pca(x, n = n, prefix = paste(prefix, "All_samples", sep = "|"), fillNA = fillNA)
    } else {
      ph <- pheno[[compare[[i]][1]]]
      sd <- setdiff(compare[[i]][-1], ph)
      if (length(sd) > 0)
        stop(sprintf("Column '%s' in 'pheno' do not contain these vaues: %s", i, paste(sd, collapse = ", ")))      
      v <- pca(x, n = n, prefix = paste(c(prefix, i), collapse = "|"), samples = which(ph %in% compare[[i]][-1]), fillNA = fillNA)
    }
    v
  })
  if (is.null(t))
    return(NULL)

  t <- t[!sapply(t, is.null)]

  list(
    samples = data.frame(do.call(cbind, lapply(t, "[[", "samples")), stringsAsFactors = FALSE, check.names = FALSE),
    features = data.frame(do.call(cbind, lapply(t, "[[", "features")), stringsAsFactors = FALSE, check.names = FALSE)
  )
}

#' Perform basic statistical analyses and preppare object to be visualized using xcmsViewer app
#' @description the PCA result is ready to be included in the prunedXcmsSet
#' @param object the prunedXcmsSet object
#' @param pheno the phenotype data.frame, could be NULL, which means pData from object will be used. 
#'   It is suggested to have a column named as "label", which contains unique value for sample labels.
#' @param median.center logical; whether to median center the data
#' @param fillNA logical; whether to fill the missing value
#' @param compare.t.test see \code{multi.t.test2}
#' @param compare.pca see \code{multi.pca}
#' @param ... other parameters passed to t.test
#' @param nf number of components to include
#' @param fx the default x-axis of feature space, e.g "ttest|A_vs_B|mean.diff", for volcano plot
#' @param fy the default y-axis of feature space, e.g "ttest|A_vs_B|log.pvalue" or "ttest|A_vs_B|log.fdr", for volcano plot
#' @param sx the default x-axis of sample space
#' @param sy the default y-axis of sample space
#' @importFrom stats p.adjust prcomp t.test
#' @export


prepViewerData <- function(
  object, pheno = NULL, median.center = FALSE, fillNA = TRUE,
  compare.t.test = NULL, compare.pca = list("_all_" = TRUE), nf = 6, ...,
  fx = "General|All|rtmed", fy = "General|All|mzmed", 
  sx = "PCA_filled|All_samples|PC1 (", sy = "PCA_filled|All_samples|PC2 ("
  ) {
  
  # phenotype data
  if (is.null(pheno))
    pd <- Biobase::pData(object@featureSet) else {
      pd <- pheno
    }
  
  expr <- Biobase::exprs(object@featureSet)
  
  if (!is.null(pd$label)) {
    pd$file <- pd$label
    rownames(pd) <- colnames(expr) <- make.names(pd$label)
  }
    
  
  # feature data
  fd <- fData(object@featureSet)
  gac <- c("ID","rtmed","mzmed","annot_ms1", "annot_ms2")
  ga <- fd[, gac]
  ea <- fd[, setdiff(colnames(fd), gac)]
  colnames(ga) <- paste("General|All", colnames(ga), sep = "|")
  colnames(ea) <- paste("General|Extended", colnames(ea), sep = "|")
  fdx <- cbind(ga, ea)
  rownames(fdx) <- rownames(fd)
  
  ######## t-test and process the data - median center / filling NAs #######
  ts <- multi.t.test2(
    x = expr, 
    compare = compare.t.test, 
    pheno = pd,
    median.center = median.center, 
    fillNA = fillNA, ...)
  
  # mat <- ts$mat
  ts <- ts$ttest
  
  ######## PCA #######
  
  pc_f <- multi.pca(expr, pheno = pd, compare = compare.pca, n = min(nrow(pd), nf), fillNA = TRUE, prefix = "PCA_filled")
  pc_nf <- multi.pca(expr, pheno = pd, compare = compare.pca, n = min(nrow(pd), nf), fillNA = FALSE, prefix = "PCA_nofill")
  
  sample_list <- list(pc_f$samples, pc_nf$samples)
  sample_list <- sample_list[!sapply(sample_list, function(x) is.null(x) || nrow(x) == 0)]

  feature_list <- list(pc_f$features, pc_nf$features)
  feature_list <- feature_list[!sapply(feature_list, function(x) is.null(x) || nrow(x) == 0)]

  pc <- list(
    samples = do.call(cbind, sample_list),
    features = do.call(cbind, feature_list)
    )

  # pc <- list(
  #   samples = cbind(pc_f$samples, pc_nf$samples),
  #   features = cbind(pc_f$features, pc_nf$features)
  #   )
  
  if (!is.null(pc$features) && nrow(pc$features) == nrow(fdx))
    fd <- cbind(fdx, pc$features) else
      fd <- fdx
  
  if (!is.null(ts))
    fd <- cbind(fd, ts[-1])
  
  colnames(pd) <- paste("General|All", colnames(pd), sep = "|")
  pd <- cbind(pd, "Stats|All|n value" = colSums(!is.na(exprs(object@featureSet))))
  if (!is.null(pc$samples) && nrow(pc$samples) == nrow(pd))
    pd <- cbind(pd, pc$samples)
  
  pData(object@featureSet) <- pd
  fData(object@featureSet) <- fd
  exprs(object@featureSet) <- expr
  
  attr(object, "fx") <- fx
  attr(object, "fy") <- fy
  attr(object, "sx") <- sx
  attr(object, "sy") <- sy
  
  object
}
