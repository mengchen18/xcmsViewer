#' Accessing elements in \code{\link[xcmsViewer]{prunedXcmsSet-class}}
#' 
#' These functions are used in the xcmsViewer shiny app. The input could be 
#' either a \code{prunedXcmsSet} object or sqlite3 database saved by 
#' \code{\link[xcmsViewer]{dumpPrunedXcmsSet}} function.
#' 
#' @section After Arguments and Value sections:
#' Despite its location, this actually comes after the Arguments and Value sections.
#' Also, don't need to use null, could annotate first function, and then
#' using function name as the groupBy name is more intuitive.
#' 
#' @param object An object of class \code{\link[xcmsViewer]{prunedXcmsSet-class}} or
#'  sqlite3 database saved by  \code{\link[xcmsViewer]{dumpPrunedXcmsSet}} 
#'  function.
#' @import RSQLite
#' @import methods
#' @name prunedXcmsSetAccesser
NULL

########### getMassTol #############
#' Access mass tolerence in PPM
#' @rdname prunedXcmsSetAccesser
#' @export
#' 
setGeneric("getMassTol", function(object) {
  standardGeneric("getMassTol")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getMassTol", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object) {
    ## need to try for back support
    v <- try(object@annot@ppmtol, silent = TRUE)
    if (inherits(v, "try-error") || is.null(v))
      return(c(MS1 = 10, MS2 = 10))
    v
  }) 

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getMassTol", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object) {
    ## need to try for back support
    v <- try(RSQLite::dbGetQuery(object, "SELECT * FROM PPMTol;"), silent = TRUE)
    if (inherits(v, "try-error") || is.null(v))
      return(c(MS1 = 10, MS2 = 10))
    structure(v$value, names = v$MSLevel)
  }) 

########### getScanIDFromFeatureID #############
#' Access feature meta data
#' @rdname prunedXcmsSetAccesser
#' @export
#' 
getScanIDFromFeatureID <- function(object, ID) { 
  fd <- getFeatureMeta(object, ID)
  pk <- getPeak(object, rowid = unlist(fd$`General|Extended|peakidx`))
  setdiff(unique(unlist(strsplit(pk$ms2Scan, ";"))), "")
}


########### getFeatureMeta #############
#' Access feature meta data
#' @rdname prunedXcmsSetAccesser
#' @export
#' 
setGeneric("getFeatureMeta", function(object, ID = NULL) {
  standardGeneric("getFeatureMeta")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getFeatureMeta", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, ID = NULL) {
    fd <- fData(object@featureSet)
    rownames(fd) <- fd[["General|All|ID"]]
    if (!is.null(ID))
      fd <- fd[fd[["General|All|ID"]] %in% ID, ]
    fd
  })

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getFeatureMeta", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object, ID = NULL) {
    if (is.null(ID)) {
      x <- RSQLite::dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_fData;")
    } else {
      x <- RSQLite::dbGetQuery(object, 'SELECT * FROM xcmsFeatureSet_fData WHERE "General|All|ID" IN (:x);', 
        param = list(x = ID))
    }
    colnames(x) <- sub("\\[", "", colnames(x))
    colnames(x) <- sub("\\]$", "", colnames(x))
    x$"General|Extended|peakidx" <- lapply(strsplit(x$"General|Extended|peakidx", ";"), as.numeric)
    rownames(x) <- x[[grep("|ID$", colnames(x))[1]]]
    x  
  }) 

########### getAnnot #############
#' Accessing annotation using monoisotopic mass
#' @param ID Metabolite feature ID
#' @rdname prunedXcmsSetAccesser
#' @export
setGeneric("getAnnot", function(object, ID) {
  standardGeneric("getAnnot")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getAnnot", 
  signature = signature(object = "prunedXcmsSet", ID = "missing"), 
  definition = function(object, ID) {
    object@annot
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getAnnot", 
  signature = signature(object = "prunedXcmsSet", ID = "character"), 
  definition = function(object, ID) {
    an <- getAnnot(object)
    .getAnnot(an, ID = ID)
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getAnnot", 
  signature = signature(object = "SQLiteConnection", ID = "missing"), 
  definition = function(object, ID) {
    RSQLite::dbGetQuery(object, "SELECT * FROM annotation;")
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getAnnot", 
  signature = signature(object = "SQLiteConnection", ID = "character"), 
  definition = function(object, ID) {
    an <- RSQLite::dbGetQuery(object, "SELECT * FROM annotation WHERE ID IN (:x);",  
                              params = list(x = ID))
    .getAnnot(an, ID = ID)
  }
)

.getAnnot <- function(an, ID) {
  an <- an[an$ID %in% ID, ]
  ms2 <- NULL
  if (nrow(an) > 0) {
    mz <- strsplit(an$ms2_mass, ";")
    intensity <- strsplit(an$ms2_intensity, ";")
    ms2 <- lapply(1:length(mz), function(i) {
      data.frame(mz = as.numeric(mz[[i]]), intensity = as.numeric(intensity[[i]]))
    })
    names(ms2) <- an$InChIKey
  }
  an$ms2_intensity <- NULL
  an$ms2_mass <- NULL
  list(ms1 = an, ms2 = ms2)
}

########### getPheno #############
#' Accessing phenotype data.frame
#' @rdname prunedXcmsSetAccesser
#' @export
setGeneric("getPheno", function(object) {
  standardGeneric("getPheno")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getPheno", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object) {
    pData(object@featureSet)
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getPheno", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object) {
    x <- RSQLite::dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_pData;")
    rownames(x) <- paste0("V", 1:nrow(x))
    x
  }
)

########### getFeatureExprs #############
#' Accessing intensity matrix/feature abundance matrix
#' @rdname prunedXcmsSetAccesser
#' @param normalized get the raw exprs or normalized exprs
#' @import Biobase
#' @export
setGeneric("getFeatureExprs", function(object, normalized = FALSE) {
  standardGeneric("getFeatureExprs")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getFeatureExprs", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, normalized = FALSE) {
    if (!normalized) {
      res <- exprs(object@featureSet)
    } else {
      res <- object@featureSet@assayData$exprs.normalized
    }
    res
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getFeatureExprs", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object, normalized = FALSE) {
    
    if (!normalized) {
      x <- RSQLite::dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_exprs;")
    } else if ("xcmsFeatureSet_exprs_normalized" %in% dbListTables(object)) {
      x <- RSQLite::dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_exprs_normalized;")
    } else {
      return(NULL)
    }
    x <- sapply(x, as.numeric)
    ids <- RSQLite::dbGetQuery(object, 'SELECT "General|All|ID" FROM xcmsFeatureSet_fData;')
    rownames(x) <- ids$'General|All|ID'
    x
  }
)
############### getPeak ##############
#' Accessing peak table
#' @param rowid the row of the peak table, set this argument if only a 
#'   subset rows of table should be returned.
#' @rdname prunedXcmsSetAccesser
#' @export
setGeneric("getPeak", function(object, rowid=NULL) {
  standardGeneric("getPeak")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getPeak", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, rowid) {
    if (is.null(rowid))
      rowid <- TRUE
    object@peak@table[rowid, ]
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getPeak", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object, rowid) {
    if (is.null(rowid)) {
      r <- RSQLite::dbGetQuery(object, "SELECT * FROM peak;")
    } else {
      r <- RSQLite::dbGetQuery(object, "SELECT * FROM peak WHERE rowid IN (:x);", params = list(x = rowid))
    }
    r
  }
)


############### getScan ##############
#' Accessing scan
#' @rdname prunedXcmsSetAccesser
#' @param scanId scan ID
#' @export
setGeneric("getScan2", function(object, scanId) {
  standardGeneric("getScan2")
})

#' @rdname prunedXcmsSetAccesser
#' @import fastmatch
#' @export
setMethod(
  "getScan2", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, scanId) {
    if (missing(scanId))
      return(object@scan@intensity)
    i <- object@scan@intensity$ID %fin% scanId
    x <- object@scan@intensity[i, ]
    a <- object@scan@meta[object@scan@meta$ID %fin% scanId, ]
    attr(x, "meta") <- a
    x
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getScan2", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object, scanId) {
    if (missing(scanId)) {
      v <- RSQLite::dbGetQuery(object, "SELECT * FROM scanIntensity;")
      return(v)
    }
    x <- RSQLite::dbGetQuery(object, "SELECT * FROM scanIntensity WHERE ID IN (:x);", params = list(x = scanId))
    a <- RSQLite::dbGetQuery(object, "SELECT * FROM scanMeta WHERE ID IN (:x);", params = list(x = scanId))
    attr(x, "meta") <- a
    x
  }
)

############### getScanMeta ##############
#' Accesssing scan meta information 
#' @rdname prunedXcmsSetAccesser
#' @export
setGeneric("getScanMeta", function(object, scanId) {
  standardGeneric("getScanMeta")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getScanMeta", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, scanId) {
    if (missing(scanId))
      return(object@scan@meta)
    i <- object@scan@meta$ID %fin% scanId
    object@scan@meta[i, ]
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getScanMeta", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object, scanId) {
    if (missing(scanId)) {
      v <- RSQLite::dbGetQuery(object, "SELECT * FROM scanMeta;")
      return(v)
    }
    RSQLite::dbGetQuery(object, "SELECT * FROM scanMeta WHERE ID IN (:x);", params = list(x = scanId))
  }
)


############### getScanUniqueID ##############
#' Accesssing unique scan IDs that have intesity after filtering
#' @rdname prunedXcmsSetAccesser
#' @export
setGeneric("getScanUniqueID", function(object) {
  standardGeneric("getScanUniqueID")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getScanUniqueID", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object) {
    unique(object@scan@intensity$ID)
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getScanUniqueID", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object) {
    RSQLite::dbGetQuery(object, "SELECT DISTINCT ID FROM scanIntensity;")[[1]]
  }
)


############### getAbundance ##############
#' Accessing feature abundance
#' @param featureId feature ID
#' @rdname prunedXcmsSetAccesser
#' @export
setGeneric("getAbundance", function(object, featureId=NULL) {
  standardGeneric("getAbundance")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getAbundance", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, featureId) {
    x <- exprs(object@featureSet)
    if (!is.null(featureId)) 
      x <- x[featureId, , drop = FALSE]
    x
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getAbundance", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object, featureId) {
    if (is.null(featureId)) {
      x <- RSQLite::dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_exprs;")
      rownames(x) <- RSQLite::dbGetQuery(object, 'SELECT "General|All|ID" FROM xcmsFeatureSet_fData;')[, "General|All|ID"] 
    } else {
      r <- RSQLite::dbGetQuery(object, 'SELECT rowid,"General|All|ID" FROM xcmsFeatureSet_fData WHERE "General|All|ID" IN (:x);', list(x = featureId))
      x <- RSQLite::dbGetQuery(object, 'SELECT * FROM xcmsFeatureSet_exprs WHERE rowid IN (:x);', list(x = r$rowid))
      rownames(x) <- r$'General|All|ID'
    }
    as.matrix(x)
  }
)

############### getEIC ##############
#' Accessing EIC
#' @param rt range of retention time
#' @param mz range of mz
#' @param file which file to extract
#' @param na.rm logical, whether remove NA values
#' @rdname prunedXcmsSetAccesser
#' @export
setGeneric("getEIC2", 
           function(object, rt = c(0, Inf), mz=c(0, Inf), 
                    file, na.rm = TRUE) {
  standardGeneric("getEIC2")
})

#' @rdname prunedXcmsSetAccesser
#' @importFrom stats na.omit
#' @export
setMethod(
  getEIC2,
  signature = signature(object = "list"), 
  definition = function(object, rt, mz, file, na.rm = TRUE) {

    mtab_subset <- object$scanMeta[object$scanMeta$msLevel == 1, ]
    mtab_subset <- mtab_subset[dplyr::between(mtab_subset$rt, rt[1], rt[2]), ]    
    if (!missing(file)) 
      mtab_subset <- mtab_subset[mtab_subset$fromFile %fin% file, ] else
        file <- unique(object$scanMeta$fromFile)
    imz <- which(object$scanIntensity$ID %fin% mtab_subset$ID & dplyr::between(object$scanIntensity$mz, mz[1], mz[2]))
    if (length(imz) == 0)
      return(
        data.frame(
          rt = numeric(0),
          intensity = numeric(0),
          file = numeric(0)
          )
        )
    
    l <- object$scanIntensity[imz, ]
    intensity <- tapply(l$intensity, l$ID, sum)
    it <- fmatch(names(intensity), object$scanMeta$ID)
    l <- data.frame(
      rt = object$scanMeta$rt[it],
      intensity = intensity,
      file = object$scanMeta$fromFile[it]
    )
    if (na.rm)
      l <- na.omit(l)
    if (length(missingfile <- setdiff(file, l$file)) > 0)
      l <- rbind(l,
                 data.frame(
                   rt = min(l$rt),
                   intensity = min(l$intensity),
                   file = missingfile
                 )
      )
    l <- l[order(l$file), ]
    l
    })

#' @rdname prunedXcmsSetAccesser
#' @importFrom stats na.omit
#' @export
setMethod(
  "getEIC2", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, rt, mz, file, na.rm = TRUE) {
    ll <- list(
      scanMeta = object@scan@meta,
      scanIntensity = object@scan@intensity
      )
    getEIC2(ll, rt = rt, mz = mz, file = file, na.rm = na.rm)
  }
)

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  "getEIC2", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object, rt, mz, file, na.rm = TRUE) {
    
    ## probably split the query to optimize the speed??
    x <- RSQLite::dbGetQuery(object, "
               SELECT m.rt,SUM(i.intensity) intensity,m.fromFile file,i.ID
               FROM scanIntensity i 
                LEFT JOIN scanMeta m
                ON i.ID=m.ID
               WHERE mz BETWEEN :u AND :v 
                AND rt BETWEEN :x AND :y 
                AND msLevel = 1
               GROUP BY i.ID
               ORDER BY fromFile;", params = list(u = mz[1], v = mz[2], x = rt[1], y = rt[2])
    )
    rownames(x) <- x$ID
    x$ID <- NULL
    if (!missing(file))
      x <- x[x$file %in% file, ]
    if (na.rm)
      x <- x[!is.na(x$intensity), ]
    if (nrow(x) == 0)
     return(NULL)
    x
  }
)

############### getEICFromFeatureID ##############
#' Accessing EIC
#' @rdname prunedXcmsSetAccesser
#' @export
setGeneric("getEICFromFeatureID", 
           function(object, ID) {
             standardGeneric("getEICFromFeatureID")
           })

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  getEICFromFeatureID,
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, ID) { 
    x <- object@EIC[[ID]]
    rownames(x) <- NULL
    x
    })

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod(
  getEICFromFeatureID,
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object, ID) { 
    RSQLite::dbGetQuery(
      object, "SELECT rt,intensity,file FROM EIC where featureID == :x;", params = list(x = ID)
      )
  })

############### getMasterPeak/not used in DB form ##############
#' Accessing master peaks
#' @rdname prunedXcmsSetAccesser
#' @export
setGeneric("getMasterPeak", function(object) {
  standardGeneric("getMasterPeak")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod("getMasterPeak", "xcmsFeatureSet", function(object) {
  assayData(object)$masterPeak
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod("getMasterPeak", "prunedXcmsSet", function(object) {
  getMasterPeak(object@featureSet)
})

############### getIonMode ##############
#' Accessing master peaks
#' @rdname prunedXcmsSetAccesser
#' @export
setGeneric("getIonMode", function(object) {
  standardGeneric("getIonMode")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod("getIonMode", "SQLiteConnection", function(object) {
  RSQLite::dbGetQuery(object, "SELECT value FROM misc WHERE param = 'IonMode';")
})

#' @rdname prunedXcmsSetAccesser
#' @export
setMethod("getIonMode", "prunedXcmsSet", function(object) {
  attr(object@annot, "IonMode")
})