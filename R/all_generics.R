################ dumpPrunedXcmsSet #############
#' Dump prunedXcmsSet to sqlite3 database
#'
#' @param object an prunedXcmsSet object
#' @param ... other arguments, such as 
#' @param db.file the directory and name of the targe database file and
#' @param overwrite overwirte the target database if existing
#' @export
#' @import RSQLite
#' @rdname prunedXcmsSet-class
#' 
setGeneric("dumpPrunedXcmsSet", function(object, ...) {
  standardGeneric("dumpPrunedXcmsSet")
})

#' @rdname prunedXcmsSet-class
setMethod(
  "dumpPrunedXcmsSet", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, db.file, overwrite = FALSE) {
    
    expr <- exprs(object@featureSet)
    cc <- pData(object@featureSet)$label
    if (is.null(cc)) {
      cc <- pData(object@featureSet)$file
      if (is.null(cc))
        cc <- rownames(cc)
    }
    colnames(expr) <- cc
    expr <- as.data.frame(expr)
    
    fd <- fData(object@featureSet)
    if (is.list(fd$peakidx))
      fd$peakidx <- sapply(fd$peakidx, paste, collapse = ";")
    
    # axis
    ax <- data.frame(
      axis = c("fx", "fy", "sx", "sy"),
      header = c(
        c(attr(object, "fx"), NA)[1],
        c(attr(object, "fy"), NA)[1],
        c(attr(object, "sx"), NA)[1],
        c(attr(object, "sy"), NA)[1]
      ),
      stringsAsFactors = FALSE
    )

    eictab <- object@EIC
    eictab <- lapply(names(eictab), function(x) {
      tab <- eictab[[x]]
      tab$featureID <- x
      tab
    })
    eictab <- do.call(rbind, eictab)

    obj <- list()
    obj$xcmsFeatureSet_fData <- fd
    obj$xcmsFeatureSet_pData <- pData(object@featureSet)
    obj$peak <- object@peak@table
    obj$scanMeta <- object@scan@meta
    obj$scanIntensity <- object@scan@intensity
    obj$annotation <- object@annot
    obj$xcmsFeatureSet_exprs <- expr
    obj$PPMTol <- attr(object@annot, "PPMTol")
    obj$defaultVis <- ax
    obj$misc <- data.frame(
      param = c("keepMS1", "IonMode"),
      value = c(object@scan@keepMS1, attr(object@annot, "IonMode"))
      )
    obj$EIC <- eictab
    # exprs.normlized
    en <- getFeatureExprs(object, normalized = TRUE)
    if (!is.null(en))
      obj$xcmsFeatureSet_exprs_normalized <- data.frame(en)
    
    ## make list as char
    obj$xcmsFeatureSet_fData$"General|Extended|peakidx" <- sapply(      
      obj$xcmsFeatureSet_fData$"General|Extended|peakidx", paste, collapse = ";"
      )
    
    mydb <- DBI::dbConnect(RSQLite::SQLite(), db.file)
    on.exit( DBI::dbDisconnect(mydb) )
    
    for (i in names(obj)) {
      cat(sprintf("Writing table %s ...\n", i))
      if (nrow(obj[[i]]) > 0) {
        DBI::dbWriteTable(mydb, name = i, value = obj[[i]], overwrite = overwrite)
      } else {
        if (overwrite) {
          rs <- RSQLite::dbSendStatement(mydb, sprintf("DROP TABLE IF EXISTS %s;", i))
          RSQLite::dbClearResult(rs)
        }
        ctc <- sprintf("CREATE TABLE %s(%s);", i, paste0("'", paste(colnames(obj[[i]]), collapse="','"), "'"))
        rs <- RSQLite::dbSendStatement(mydb, ctc, overwrite = overwrite)
        RSQLite::dbClearResult(rs)
      }
    }
    
    cat("Creating indices ...\n")
    indices <- c(
      "CREATE UNIQUE INDEX index_fData_ID ON xcmsFeatureSet_fData ('General|All|ID');",
      "CREATE UNIQUE INDEX index_scanMeta_ID ON scanMeta (ID);",
      "CREATE UNIQUE INDEX index_scanMeta_rt_ID_file_msLevel ON scanMeta (rt, ID, fromFile, msLevel);",
      "CREATE INDEX index_scanIntensity_ID ON scanIntensity (ID);",
      "CREATE INDEX index_scanIntensity_mz ON scanIntensity (mz);",
      "CREATE INDEX index_eic_fid ON EIC (featureID);"
    )
    for (stat in indices) {
      res <- RSQLite::dbSendStatement(mydb, stat)
      RSQLite::dbClearResult(res)
    }
  })

######################## annotateMetabolite ########################
#' Annotate prunedXcmsSet
#'
#' @param mode Ionization mode. pos or neg
#' @param ref reference MS2 spectra used annotate the experimental MS2 spectra
#' @param ppmtol the mass tolerence given by parts per million (PPM)
#' @param ips ips
#' @param fun_parallel parallel function, often either mclapply or biocparallel
#' @export
#' @rdname prunedXcmsSet-class
#' 
setGeneric("annotateMetabolite", function(object, ...) {
  standardGeneric("annotateMetabolite")
})

#' @rdname prunedXcmsSet-class
setMethod(
  "annotateMetabolite", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(
    object, mode = c('pos', "neg")[1], ref, ppmtol = 25, ips = 0.75, fun_parallel= parallel::mclapply, ...) {
    
    
    x <- object    
    fd <- fData(x@featureSet)    
    fd$annot_ms2 <- NA
    ptol <- data.frame(
      MSLevel = c("MS1", "MS2"),
      value = c(ppmtol, ppmtol),
      stringsAsFactors = FALSE)
    
    an1 <- annotateMS1(
      x, mode = mode, ref = ref, ppmtol = ppmtol, ips = ips, fun_parallel = fun_parallel, ...
    )    
    an2 <- annotateMS2(x = x, ms1Annot = an1, ppmtol = ppmtol)    
    an <- scoreAnnot(x = x, an2 = an2)
    an <- an[order(an$Score, decreasing = TRUE),]    
    an <- an[.xcmsViewerInternalObjects()$xcmsAnnot_column]        
    r <- categorizeAnnotation(fd = fd, an = an, ppmtol = ppmtol)
    an <- r$an
    attr(an, "PPMTol") <- ptol    
    attr(an, "IonMode") <- mode    

    fd <- r$fd
    ms1 <- sapply(unique(an$ID), function(x) {
      nam <- an$cpdName[an$ID == x]
      nam <- nam[1:min(3, length(nam))]
      substr(paste(nam, collapse = ";"), 1, 60)
    })    
    ms1 <- ms1[rownames(fd)]    
    ms1[is.na(ms1)] <- ""    
    fd$annot_ms1 <- ms1    
    
    fData(x@featureSet) <- fd    
    x@annot <- an    
    x
  })


.validRestrictedDataFrame <- function(df, name, class, str, contain = FALSE) {
  if (ncol(df) != 0)  {
    if (contain) {
      i <- setdiff(name, colnames(df))
      if (length(i) > 0) 
        return(sprintf("Missing required columns in %s: %s", str, paste(i, collapse = ", ")))
      df <- df[, name]
    } else {
      if (!identical(colnames(df), name))
        return(sprintf("Problem in %s column name!", str))
    }
    if (!missing(class)) {
      ii <- sapply(seq_along(name), function(x) {
        inherits(df[[x]], c(class[x], "AsIs"))
        })
      # if (!all( sapply(df, class) == class))
      if (!all( ii ))      
        return(sprintf("Problem in %s column class!", str))
    }
  }
}

######################## exportTables ########################
#' Export metabolite and pheno table
#'
#' @param file xlsx file
#' @export
#' @import openxlsx
#' @importFrom Biobase exprs pData fData
#' @rdname prunedXcmsSet-class

setGeneric("exportTables", function(object, file) {
  standardGeneric("exportTables")
})

#' @rdname prunedXcmsSet-class
setMethod(
  "exportTables", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, file) { 
    pd <- Biobase::pData(object@featureSet)    
    fd <- Biobase::fData(object@featureSet)
    fd$peakidx <- sapply(fd$peakidx, paste, collapse = ",")
    expr <- Biobase::exprs(object@featureSet)
    if ( "label" %in% colnames(pd) ) 
      coln <- pd$label else
        coln <- pd[, 1]
    colnames(expr) <- coln
    .writeMetaboliteTables(pd, fd, expr, file)
  })

#' @rdname prunedXcmsSet-class
setMethod(
  "exportTables", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object, file) { 
    pd <- dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_pData")
    fd <- dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_fData")
    expr <- dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_exprs_normalized")
    if ( "label" %in% colnames(pd) ) 
      coln <- pd$label else
        coln <- pd[, 1]
    colnames(expr) <- coln
    .writeMetaboliteTables(pd, fd, expr, file)
  })

.writeMetaboliteTables <- function(pd, fd, expr, file) {
  wb <- createWorkbook("BayBioMS")
  addWorksheet(wb, sheetName = "SampleInfo")
  addWorksheet(wb, sheetName = "MetaboliteInfo")
  addWorksheet(wb, sheetName = "MetaboliteIntensity")
  writeData(wb, sheet = "SampleInfo", x = pd)
  writeData(wb, sheet = "MetaboliteInfo", x = fd)
  writeData(wb, sheet = "MetaboliteIntensity", x = cbind(ID = fd$ID, expr))
  saveWorkbook(wb, file = file, overwrite = TRUE)
}


