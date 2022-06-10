
categorizeAnnotation <- function(fd, an, ppmtol) {
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ function start ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  removePrecusorPeak <- function(x, rm.mass, ppmtol = 30) {
    rm.mass <- rm.mass - rm.mass*ppmtol/1e6
    x[which(x$mz < rm.mass), ]
  }
  
  ic_ms2 <- grep("ms2spectrum$", colnames(fd))
  
  # ~~~~~~ iterate from here ~~~~~~~~
  lvs <- vapply(seq_len(nrow(an)), function(i) {
    an1 <- an[i, ]
    ifd <- which(fd$ID == an1$ID)
    
    str_ms2_annot <- an1$ms2_mass
    str_ms2_measured <- fd[ifd, ic_ms2]
    rt_annot <- an1$RT
    rt_measured <- fd$rtmed[ifd]
    
    spec_annot <- data.frame(
      mz = as.numeric(strsplit(an1$ms2_mass, split = ";")[[1]]),
      intensity = as.numeric(strsplit(an1$ms2_intensity, split = ";")[[1]]),
      stringsAsFactors = FALSE
    )
    spec_measured <- str2spectra( fd$ms2spectrum[ifd] )
    
    
    spec_annot <- removePrecusorPeak(spec_annot, rm.mass = an1$MassWithAdduct)
    spec_measured <- removePrecusorPeak(spec_measured, rm.mass = fd$mzmed[ifd])
    
    if (is.null(spec_measured) || nrow(spec_annot) < 1 || nrow(spec_measured) < 1) {
      if (is.na(rt_annot) || nchar(rt_annot) < 1) {
        lev <- 1
      } else {
        rtd <- abs(rt_measured - rt_annot)
        if (rtd < 30)
          lev <- 3 else if (rtd > 90)
            lev <- 0 else
              lev  <- 1
      }
    } else {
      pmap <- .prep_mirrorPlot(spec_measured, spec_annot, ppmtol = ppmtol)
      int <- abs(pmap$tab$intensity)
      ep <- sum(int[!is.na(pmap$tab$mapGroup)])/sum(int)
      if (ep < 0.05)  {
        lev <- 0
      } else if (ep < 0.75) {
        lev <- 2
      } else 
        lev <- 4
    }
    as.integer(lev)
  }, FUN.VALUE = integer(1))
  an <- data.frame( category = c(lvs), an )
  an <- an[order(an$category, decreasing = TRUE), ]
  
  # update fd
  flv <- tapply(an$category, an$ID, max)
  flv <- c(flv[fd$ID])
  flv[is.na(flv)] <- "NI"
  fd$"annot_ms2" <- flv
  
  list(an = an, fd = fd)
}
