#' Query possible metabolites according to given m/z
#' @param m the mass of the metabolites
#' @param tolppm the mass tolerence given by parts per million (PPM)
#' @param refTab the reference metabolites table, it must contain at list the following columns:
#'   name - the name of metabolites
#'   monoisotopic_molecular_weight - the monoisotopic molecular weight of metabolites
#' @param addTable the adduct table, it must contain at least the following columns:
#'   adduct - the name of the adduct, such as "[M+H]+"
#'   massdiff - the mass different of an adduct
#'   nmol - the mol of metabolite, for example, [2M-H] will be 2
#'   ips - Four values are possible (0.25,0.5,0.75,1) depending on the likelihood of the rule 
#' @param ID an optional ID column added at the end of the matched data.frame
#' @export

massQuery <- function(m, tolppm, refTab, addTable, ID=NULL) {
  
  if (is.null(refTab$monoisotopic_molecular_weight))
    stop("The 'monoisotopic_molecular_weight' column doesn't exist in the refTab!")
  
  m_mass <- (m - addTable$massdiff)/addTable$nmol
  dd <- sapply(m_mass, function(v) {
    v - refTab$monoisotopic_molecular_weight
  })
  
  i <- which(abs(dd) < tolppm*m/1e6, arr.ind = TRUE)
  if (nrow(i) == 0)
    return(NULL)
  at <- cbind(refTab[i[, 1], ], addTable[i[, 2], ])
  at$massWithAdduct <- at$monoisotopic_molecular_weight + at$massdiff
  at$massQueried <- m
  at$deltaPPM <- abs(at$massQueried - at$massWithAdduct)/at$massQueried * 1e6
  at$score <- (1-at$deltaPPM/tolppm) * at$ips
  at <- at[order(at$score, decreasing = TRUE), ]
  at$ID <- ID
  at
}
