dir_fun <- "/media/share_baybioms/Projects/008_Bioinformatics/B012_MAIT_devel/xcmsViewerData/R"
dir_ms_dial <- "/media/share_baybioms/Lab_organisation/Massenspektrometer/Useful-Tools/MSP-Library_from-MSDIAL/"
dir_rds <- "/media/share_baybioms/Projects/008_Bioinformatics/B012_MAIT_devel/xcmsViewerData/inst/extdata/"
hmdb_ms1_x <- "/media/share_baybioms/Projects/008_Bioinformatics/B012_MAIT_devel/Databases/HMDB/Rds/hmdb_metabolites.rds"
hmdb_ms2_exp <- "/media/share_baybioms/Projects/008_Bioinformatics/B012_MAIT_devel/Databases/HMDB/MS2peaks/hmdb_experimental_msms_spectra/"
hmdb_ms2_pred <- "/media/share_baybioms/Projects/008_Bioinformatics/B012_MAIT_devel/Databases/HMDB/MS2peaks/hmdb_predicted_msms_spectra/"


# =================================================================================
#
#                ms dial
#
# =================================================================================

fls <- list.files(dir_ms_dial, pattern = ".msp$")
for (i in fls) {
  print(i)
  dbn <- sub(".msp$", "", i)
  dbn <- sub("^MSMS-", "", dbn)
  db <- prep_annotation_msdial(file.path(dir_ms_dial, i), db = dbn, predicted = "n/a", cores = 60)
  saveRDS(db$meta_ms1, file.path(dir_rds, "an_ms1", paste0("MSDIAL_", dbn, ".rds")))
  saveRDS(list(meta = db$meta_ms2, peakList = db$peaks), 
          file.path(dir_rds, "an_ms2", paste0("MSDIAL_", dbn, ".rds")))
}



# =================================================================================
#
#                hmdb
#
# =================================================================================

hmdb_ms1 <- prep_annotation_hmdb(hmdb_ms1_x)
saveRDS(hmdb_ms1, file = file.path(dir_rds, "an_ms1", "HMDB.rds"))
gc()


library(stringr)
library(parallel)
library(flatxml)
library(parallel)
library(BiocParallel)

.laf <- list.files(dir_fun, full.names = TRUE)
for (.i in .laf)
  source(.i)

hmdb_ms2_exp <- prep_annotation_hmdb_msms(xmlFolder = hmdb_ms2_exp)
saveRDS(hmdb_ms2_exp, file = file.path(dir_rds, "an_ms2", "HMDB_experimental.rds"))
gc()

hmdb_ms2_pred <- prep_annotation_hmdb_msms(xmlFolder = hmdb_ms2_pred)
saveRDS(hmdb_ms2_pred, file = file.path(dir_rds, "an_ms2", "HMDB_predicted.rds"))
gc()

