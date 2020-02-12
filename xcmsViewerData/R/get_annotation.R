#' Get annotation of masses
#' @param what the annotation to get, the options are:
#'    HMDB
#'    MSDIAL_BioMSMS-Neg-PlaSMA
#'    MSDIAL_BioMSMS-Pos-PlaSMA
#'    MSDIAL_KI-GIAR_zic-HILIC_Pos_v0.90
#'    MSDIAL_Neg-CASMI2016
#'    MSDIAL_Neg-FiehnHILIC
#'    MSDIAL_Neg-GNPS
#'    MSDIAL_Neg-MassBankEU
#'    MSDIAL_Neg-MassBank
#'    MSDIAL_Neg-MetaboBASE
#'    MSDIAL_Neg-PlaSMA
#'    MSDIAL_Neg-Respect
#'    MSDIAL_Neg-RikenOxPLs
#'    MSDIAL_Neg-Vaniya-Fiehn_Natural_Products_Library_20200109
#'    MSDIAL_Pos-CASMI2016
#'    MSDIAL_Pos-FiehnHILIC
#'    MSDIAL_Pos-GNPS
#'    MSDIAL_Pos-MassBankEU
#'    MSDIAL_Pos-MassBank
#'    MSDIAL_Pos-MetaboBASE
#'    MSDIAL_Pos-Pathogen_Box_20200109
#'    MSDIAL_Pos-PlaSMA
#'    MSDIAL_Pos-Respect
#'    MSDIAL_Pos-Vaniya-Fiehn_Natural_Products_Library_20200109
#'    MSDIAL_Public-Neg-VS14
#'    MSDIAL_Public-Pos-VS14
#' @export
#' 

get_annotation_mass <- function(what) {
  
  opts <- c(
    "HMDB", 
    "MSDIAL_BioMSMS-Neg-PlaSMA", 
    "MSDIAL_BioMSMS-Pos-PlaSMA", 
    "MSDIAL_KI-GIAR_zic-HILIC_Pos_v0.90", 
    "MSDIAL_Neg-CASMI2016", 
    "MSDIAL_Neg-FiehnHILIC", 
    "MSDIAL_Neg-GNPS", 
    "MSDIAL_Neg-MassBankEU", 
    "MSDIAL_Neg-MassBank", 
    "MSDIAL_Neg-MetaboBASE", 
    "MSDIAL_Neg-PlaSMA", 
    "MSDIAL_Neg-Respect", 
    "MSDIAL_Neg-RikenOxPLs", 
    "MSDIAL_Neg-Vaniya-Fiehn_Natural_Products_Library_20200109", 
    "MSDIAL_Pos-CASMI2016", 
    "MSDIAL_Pos-FiehnHILIC", 
    "MSDIAL_Pos-GNPS", 
    "MSDIAL_Pos-MassBankEU", 
    "MSDIAL_Pos-MassBank", 
    "MSDIAL_Pos-MetaboBASE", 
    "MSDIAL_Pos-Pathogen_Box_20200109", 
    "MSDIAL_Pos-PlaSMA", 
    "MSDIAL_Pos-Respect", 
    "MSDIAL_Pos-Vaniya-Fiehn_Natural_Products_Library_20200109", 
    "MSDIAL_Public-Neg-VS14", 
    "MSDIAL_Public-Pos-VS14")
  if (any(!what %in% opts))
    "Unknown resource!"
  
  dir <- system.file(package = "xcmsViewerData", "extdata", "an_ms1")
  dl <- lapply(what, function(w) {
    readRDS(file.path(dir, paste0(w, ".rds")))
  })
  do.call(rbind, dl)
}


#' Get annotation of fragments (MS2)
#' @param what the annotation to get, the options are:
#'    HMDB_experimental
#'    HMDB_predicted
#'    MSDIAL_BioMSMS-Neg-PlaSMA
#'    MSDIAL_BioMSMS-Pos-PlaSMA
#'    MSDIAL_KI-GIAR_zic-HILIC_Pos_v0.90
#'    MSDIAL_Neg-CASMI2016
#'    MSDIAL_Neg-FiehnHILIC
#'    MSDIAL_Neg-GNPS
#'    MSDIAL_Neg-MassBankEU
#'    MSDIAL_Neg-MassBank
#'    MSDIAL_Neg-MetaboBASE
#'    MSDIAL_Neg-PlaSMA
#'    MSDIAL_Neg-Respect
#'    MSDIAL_Neg-RikenOxPLs
#'    MSDIAL_Neg-Vaniya-Fiehn_Natural_Products_Library_20200109
#'    MSDIAL_Pos-CASMI2016
#'    MSDIAL_Pos-FiehnHILIC
#'    MSDIAL_Pos-GNPS
#'    MSDIAL_Pos-MassBankEU
#'    MSDIAL_Pos-MassBank
#'    MSDIAL_Pos-MetaboBASE
#'    MSDIAL_Pos-Pathogen_Box_20200109
#'    MSDIAL_Pos-PlaSMA
#'    MSDIAL_Pos-Respect
#'    MSDIAL_Pos-Vaniya-Fiehn_Natural_Products_Library_20200109
#'    MSDIAL_Public-Neg-VS14
#'    MSDIAL_Public-Pos-VS14
#' @export

get_annotation_fragment <- function(what) {
  
  opts <- c("HMDB_experimental", 
            "HMDB_predicted", 
            "MSDIAL_BioMSMS-Neg-PlaSMA", 
            "MSDIAL_BioMSMS-Pos-PlaSMA", 
            "MSDIAL_KI-GIAR_zic-HILIC_Pos_v0.90", 
            "MSDIAL_Neg-CASMI2016", 
            "MSDIAL_Neg-FiehnHILIC", 
            "MSDIAL_Neg-GNPS", 
            "MSDIAL_Neg-MassBankEU", 
            "MSDIAL_Neg-MassBank", 
            "MSDIAL_Neg-MetaboBASE", 
            "MSDIAL_Neg-PlaSMA", 
            "MSDIAL_Neg-Respect", 
            "MSDIAL_Neg-RikenOxPLs", 
            "MSDIAL_Neg-Vaniya-Fiehn_Natural_Products_Library_20200109", 
            "MSDIAL_Pos-CASMI2016", 
            "MSDIAL_Pos-FiehnHILIC", 
            "MSDIAL_Pos-GNPS", 
            "MSDIAL_Pos-MassBankEU", 
            "MSDIAL_Pos-MassBank", 
            "MSDIAL_Pos-MetaboBASE", 
            "MSDIAL_Pos-Pathogen_Box_20200109", 
            "MSDIAL_Pos-PlaSMA", 
            "MSDIAL_Pos-Respect", 
            "MSDIAL_Pos-Vaniya-Fiehn_Natural_Products_Library_20200109", 
            "MSDIAL_Public-Neg-VS14", 
            "MSDIAL_Public-Pos-VS14")
  
  if (any(!what %in% opts))
    "Unknown resource!"
  
  dir <- system.file(package = "xcmsViewerData", "extdata", "an_ms2")
  dl <- lapply(what, function(w) {
    readRDS(file.path(dir, paste0(w, ".rds")))
  })
  combine_ms2(dl)
}