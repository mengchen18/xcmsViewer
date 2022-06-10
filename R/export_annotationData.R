#' Database used to annotate identified features in untargeted metabolomics experiment
#'
#' The dataset is compiled from two resources: 1) The experimental MS-MS sepctra 
#'   from HMDB and 2) all public MS/MS from MSDIAL. 
#'
#' \itemize{
#'   \item InChIKey. InChIKey of metabolite
#'   \item CID. Pubchem compound ID.
#'   \item cpdName. Compound name.
#'   \item formula. Compound formula. 
#'   \item monoMass. Monoisotopic mass of metabolite.
#'   \item POS_mass. Peak mass of consensus MS/MS spectra in positive ionization mode.
#'   \item POS_intensity. Peak intensity of consensus MS/MS spectra in positive ionization mode.
#'   \item POS_purity. Purity of consensus MS/MS spectra in positive ionization mode. 
#'     It should be a value between 0-1, 1 means the individual MS/MS spectrum has good agreement. 
#'   \item POS_sourceId. The IDs of combined MS/MS spectra combined to get 
#'     the consensus MS/MS spectra, in positive mode. 
#'   \item NEG_mass. See POS_mass, for negative ionization mode. 
#'   \item NEG_intensity. See POS_intensity, for negative ionization mode. 
#'   \item NEG_purity. See POS_purity, for negative ionization mode. 
#'   \item NEG_sourceId. See POS_sourceId, for negative ionization mode. 
#'   \item smiles. SMILES of metabolite. 
#'   \item RT. Retention time. For public database, this column is set to NA. User can add in-house 
#'     annotations and give RT so the RT will also be considered in metabolite annotation. 
#' }
#'
#' @docType data
#' @keywords datasets
#' @name hmdb_msdial
#' @usage data(hmdb_msdial)
#' @format A data frame with 120000+ rows and 14 variables
#' @source <https://hmdb.ca/downloads>
#' @source <http://prime.psc.riken.jp/compms/msdial/main.html#MSP>
NULL
