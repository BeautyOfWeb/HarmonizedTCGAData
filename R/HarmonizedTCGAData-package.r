#' HarmonizedTCGAData.
#'
#' @name HarmonizedTCGAData
#' @docType package
NULL


#' Wall
#'
#' Wall contains a list of precomputed affinity (similarity) matrices of 2582 patients.
#' These matrices were derived from 10382 gene expression, miRNA expression and
#' DNA methylation data files downloaded from GDC data portal
#' The file UUIDs can be found in inst/extdata/fileUUIDs.csv
#' Using these file UUIDs, users can download the original data from
#' https://portal.gdc.cancer.gov/repository
#' `Wall` is a complex list and contains lists inside list.
#' Precisely, Wall is a list (five cancer types) of list (six feature normalization types:
#' raw.all, raw.sel, log.all, log.sel, vst.sel, normalized) of list (three feature spaces
#' or views: fpkm, mirna, and methy450) of matrices.
#' (So Wall contains 90 matrices in total)
#' The rownames of each matrix is the case_id (i.e., patient id),
#' and the column names of each matrix is the aliquot IDs (i.e., TCGA barcode,
#' which contains the case_id as prefix).
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' myfiles <- query(eh, "HarmonizedTCGAData")
#' Wall <- myfiles[[1]]
#' # Wall <- myfiles[['EH1014']]
#' names(Wall)
#' names(Wall[[1]])
#' names(Wall[[1]][[1]])
#' dim(Wall[[1]][[1]][[1]])
"Wall"


#' project_ids
#'
#' A named character vector: mapping case_id (i.e., patient ID) to
#' the TCGA project ID they belong to
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' myfiles <- query(eh, "HarmonizedTCGAData")
#' project_ids <- myfiles[[2]]
#' # project_ids <- myfiles[['EH1015']]
#' head(project_ids)
"project_ids"



#' surv.plot
#'
#' Patient survival information (overall survival plot data) were downloaded from
#' https://portal.gdc.cancer.gov/exploration?searchTableTab=genes
#' For detailed information: see section "Survival analysis" in
#' https://docs.gdc.cancer.gov/Data_Portal/PDF/Data_Portal_UG.pdf
#' @format A data frame with four variables: \code{survivalEstimate}, \code{id},
#'   \code{censored}, and \code{time}
#' @import ExperimentHub
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' myfiles <- query(eh, "HarmonizedTCGAData")
#' surv.plot <- myfiles[[3]]
#' # surv.plot <- myfiles[['EH1016']]
#' head(surv.plot)
"surv.plot"
