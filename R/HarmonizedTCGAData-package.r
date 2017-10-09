#' HarmonizedTCGAData.
#'
#' @name HarmonizedTCGAData
#' @docType package
NULL


#' surv.plot
#'
#' Patient survival information downloaded from https://portal.gdc.cancer.gov
#' For detailed information: see section "Survival analysis" in https://docs.gdc.cancer.gov/Data_Portal/PDF/Data_Portal_UG.pdf
#' @format A data frame with four variables: \code{survivalEstimate}, \code{id},
#'   \code{censored}, and \code{time}
"surv.plot"

#' project_ids
#' A named character vector: mapping sumbitter_id (can be seen as patient ID) to its TCGA project ID
"project_ids"

#' Wall
#' Wall is a complex list and contains lists inside list.
#' Precisely, Wall is a list (five cancer types) of list (six feature normalization types: raw.all, raw.sel, log.all, log.sel, vst.sel, normalized) of list (three feature spaces or views: fpkm, mirna, and methy450) of matrices.
#' The rownames of each matrix is the submitter_id (can be seen as a patient id),
#' and the column names of each matrix is the aliquot ID (which contains the submitter_id as prefix).
#' Based on these aliquot IDs, users can download original data from https://portal.gdc.cancer.gov/repository.
"Wall"
