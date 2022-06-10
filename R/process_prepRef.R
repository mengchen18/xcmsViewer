#' preprare reference 
#' @param ref annotation data.frame
#' @param mode mode, either pos or neg
#' @param primaryInChI the InChIkey for primary annotation, not used yet, so 
#'   users need to edit the ref before passed to the runPrunedXcmsSet beforehand.
#' 
prepRef <- function(ref, mode, primaryInChI = NULL) {
  
  mode <- tolower(mode)
  mode <- match.arg(mode, c("positive", "negative"))
  if (mode == "positive") {
    ref$ms2_mass <- ref$POS_mass
    ref$ms2_intensity <- ref$POS_intensity
    ref$ms2_purity <- ref$POS_purity
    ref$ms2_sourceId <- ref$POS_sourceId
  } else if (mode == "negative") {
    mod <- "neg"
    ref$ms2_mass <- ref$NEG_mass
    ref$ms2_intensity <- ref$NEG_intensity
    ref$ms2_purity <- ref$NEG_purity
    ref$ms2_sourceId <- ref$NEG_sourceId
  } else {
    stop("Unknown mode!")
  }
  
  if (!"RT" %in% colnames(ref))
    ref$RT <- NA
  if (! "primary" %in% colnames(ref))
    ref$primary <- FALSE
  if (!is.null(primaryInChI))
    ref$primary <- ref$InChIKey %in% primaryInChI
  
  ref
}