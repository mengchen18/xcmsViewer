#' xcms annotation of MS1
#' @param x an object of prunedXcmsSet
#' @param mode Ionization mode. pos or neg
#' @param ref reference MS2 spectra used annotate the experimental MS2 spectra
#' @param ppmtol the mass tolerence given by parts per million (PPM)
#' @param fun_parallel the parallel function, could be \code{mclapply} or \code{bplapply}
#' @param ips The likelihood of the rule ("ips" column) in adduct table in MAIT package. 
#' @param ... the parameters passed to parallel function
#' @import parallel
#' @import BiocParallel
#' @importFrom MAIT Biotransformations
#' @importFrom utils data
#' 
annotateMS1 <- function(
  x, mode = c("pos", "neg")[1], ref, ppmtol = 10, ips = 0.75, fun_parallel = parallel::mclapply, ...
) {
  features <- fData(x@featureSet)
  maitEnv <- environment()
  data("MAITtables", package = "MAIT", envir = maitEnv)
  if (mode == "pos") {
    at <- maitEnv$posAdducts    
  } else if (mode == "neg") {
    at <- maitEnv$negAdducts
  } else 
      stop("'mode' should be either pos or neg!")
  at <- data.frame(
    Adduct = as.character(at$name),
    MassDiff = at$massdiff,
    Nmol = as.numeric(at$nmol),
    IPS = as.numeric(at$ips),
    stringsAsFactors = FALSE)
  at <- at[which(at$IPS >= ips), ]
  
  ll <- fun_parallel(1:nrow(features), function(i) {
    massQuery(m = features$mzmed[i], ppmtol = ppmtol, refTab = ref, 
              addTable = at, ID = features$ID[i]) 
  }, ...)
  v <- do.call(rbind, ll)
  attr(v, "ppmtol") <- ppmtol
  v
}


#' Query possible metabolites according to given m/z
#' @param m the mass of the metabolites
#' @param ppmtol the mass tolerence given by parts per million (PPM)
#' @param refTab the reference metabolites table, it must contain at list the following columns:
#'   name - the name of metabolites
#'   monoisotopic_molecular_weight - the monoisotopic molecular weight of metabolites
#' @param addTable the adduct table, it must contain at least the following columns:
#'   adduct - the name of the adduct, such as "[M+H]+"
#'   massdiff - the mass different of an adduct
#'   nmol - the mol of metabolite, for example, [2M-H] will be 2
#'   ips - Four values are possible (0.25,0.5,0.75,1) depending on the likelihood of the rule 
#' @param ID an optional ID column added at the end of the matched data.frame
#' @examples
#' # library(MAIT)
#' # library(xcmsViewerData)
#' # data(MAITtables)
#' # m1 <- get_annotation_mass("HMDB")
#' # q <- massQuery(168.07864 * 2 + 1.007276, ppmtol=10, refTab=m1, addTable=posAdducts, ID=NULL)
#' @export

massQuery <- function(m, ppmtol, refTab, addTable, ID=NULL) {
  
  if (is.null(refTab$monoMass))
    stop("The 'monoMass' column doesn't exist in the refTab!")
  
  m_mass <- (m - addTable$MassDiff)/addTable$Nmol
  dd <- sapply(m_mass, function(v) {
    v - refTab$monoMass
  })
  
  i <- which(abs(dd) < ppmtol*m/1e6, arr.ind = TRUE)
  if (nrow(i) == 0)
    return(NULL)
  
  at <- data.frame(ID = rep(ID, nrow(i)), stringsAsFactors = FALSE)
  at <- cbind(at, refTab[i[, 1], c("InChIKey","CID","cpdName","formula", "monoMass")], addTable[i[, 2], ])
  at$MassWithAdduct <- at$monoMass + at$MassDiff
  at$MassQueried <- m
  at$DeltaPPM <- abs(apply(i, 1, function(x) dd[x[1], x[2]])/m*1e6)
  at <- cbind(at, refTab[i[, 1], -(1:4)])
  at
}
