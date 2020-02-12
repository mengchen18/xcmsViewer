.ms_meta_header <- function() {
  c("name", "chemical_formula", "average_molecular_weight", 
  "monisotopic_molecular_weight", "cas_registry_number", "accession")
}

.msms_meta_header <- function() {
  c("id","notes","sample-concentration","solvent","sample-mass",
    "sample-assessment","spectra-assessment","sample-source","collection-date","instrument-type",
    "peak-counter","created-at","updated-at","mono-mass","energy-field",
    "collision-energy-level","collision-energy-voltage","ionization-mode",
    "sample-concentration-units","sample-mass-units",
    "predicted","structure-id","splash-key","database-id")
}

#' id creater for MS2
#' @param num the number 
#' @param db the db prefix
#' @importFrom stringr str_pad
 
.msms_id_creater <- function(num, db) {
  paste0(db, "_MSMS", stringr::str_pad(num, width = 9, pad = "0"))
}

#' get value from a named vector, if not, return NA
#' @param x a named vector
#' @param val value to get
#' @importFrom fastmatch fmatch

.getVal <- function(x, val) {
  i <- fmatch(val, names(x))
  if (is.na(i))
    return(i)
  x[[i]]
}



#' Convert precusor MZ to monosiotopic mass
#' @param type type of adduct, such as "[M+Cl]-"
#' @param precMass the precursor mass
#' @import MAIT
#' @examples
#' # q <- .toMonoisoMass(type = c("[M-H]-", "[M-H]-"), precMass = c(194.91767134783, 174.97229374783))
#' # c(2,3,5-Trichlorphenol, 2,4-Dichloro-6-Methylphenol)
#' 
.toMonoisoMass <- function(type, precMass) {
  
  trans <- c("M+Cl" = "[M+Cl]-",
             "M+K" =  "[M+K]+",
             "M+NH4" = "[M+NH4]+",
             "M+Na" = "[M+Na]+",
             "M+H" = "[M+H]+",
             "M-H" = "[M-H]-" )
  
  ir <- fastmatch::fmatch(type, names(trans))
  ir <- trans[ir]
  ina <- !is.na(ir)
  type[ina] <- ir[ina]
  
  tabmaitenv <- new.env()
  load(system.file(package = "MAIT", "data", "MAITtables.RData"), envir = tabmaitenv)
  
  negAdd <- tabmaitenv$negAdducts
  negAdd <- negAdd[negAdd$nmol == 1, ]
  negAdd <- structure(negAdd$massdiff, names = as.character(negAdd$name))
  
  posAdd <- tabmaitenv$posAdducts
  posAdd <- posAdd[posAdd$nmol == 1, ]
  posAdd <- structure(posAdd$massdiff, names = as.character(posAdd$name))
  add <- c(negAdd, posAdd)
  
  precMass - add[type]
}





