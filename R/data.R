#' Example expression dataset for REBOOT
#'
#' @description
#' A toy expression dataset included in the REBOOT package intended for demonstration, testing, and vignette examples.
#'
#' @details
#' The dataset contains expression data along with survival information.
#' Each row represents a sample (patient), and each column represents a feature (gene/transcript)
#'
#' The first two columns correspond to survival data:
#' \itemize{
#'   \item \code{OS}: Overall survival status (0 = censored, 1 = event)
#'   \item \code{OS.time}: Overall survival time
#' }
#'
#' The remaining columns correspond to expression values.
#'
#' Row names represent unique sample identifiers (e.g., patient IDs).
#'
#' @format
#' A \code{data.frame} with:
#' \itemize{
#'   \item Rows: samples (patients)
#'   \item Columns: survival variables + expression features
#' }
#'
#' @source
#' Derived from internal dataset (FREDY REF here), reduced for package examples.
#'
#' @usage data(toy_expression)
#'
#' @examples
#' data(toy_expression)
#' dim(toy_expression)
#' head(toy_expression[, 1:5])
#'
#' @keywords datasets
"toy_expression"

#' Example clinical dataset for REBOOT
#'
#' @description
#' A toy clinical dataset of Brain Lower Grade Glioma (TCGA-LGG) patients included in the REBOOT package
#' intended to complement the expression dataset with additional sample annotations.
#'
#' @details
#' Each row represents a sample (patient), matching the expression dataset.
#'
#' The dataset includes the following categorical binary clinical variables:
#' \itemize{
#'   \item Age group (e.g., Low/High)
#'   \item Gender
#'   \item Race (e.g., White/Not white)
#'   \item Histological type (e.g., Astrocytoma/Oligodendroglioma)
#'   \item Tumor grade (e.g., G2/G3)
#' }
#'
#' Row names correspond to TCGA patient IDs and are aligned with \code{toy_expression}.
#'
#' @format
#' A \code{data.frame} with:
#' \itemize{
#'   \item Rows: samples (patients)
#'   \item Columns: clinical variables
#' }
#'
#' @source
#' Derived from internal dataset (FREDY REF here), reduced for package examples.
#'
#' @usage data(toy_clinics)
#'
#' @examples
#' data(toy_clinics)
#' dim(toy_clinics)
#' head(toy_clinics)
#'
#' @keywords datasets
"toy_clinics"
