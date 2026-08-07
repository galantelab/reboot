###### CONSTRUCTORS ######

#' Internal constructor for REBOOT signature object
#'
#' @description Creates a standardized `reboot_signature` S3 object from a signature table. This function ensures that
#' the input signature is valid, properly formatted, and consistent across REBOOT modules.
#'
#' @details This constructor is used internally by REBOOT functions (e.g., \code{rebootRegression()}
#' and \code{rebootSurvival()}) to guarantee object integrity and standardization across REBOOT modules.
#'
#' It performs the following steps:
#' \itemize{
#'   \item Validates that the input is a data.frame
#'   \item Ensures required columns (\code{feature}, \code{coefficient}) are present
#'   \item Removes rows with missing coefficients
#'   \item Sorts features by absolute coefficient value (descending)
#' }
#'
#' Additional components such as plots, parameters, metadata, and function call can be optionally attached to the object.
#'
#' @param signature A data.frame containing at least two columns:
#'   \itemize{
#'     \item \code{feature}: Feature names (e.g., genes or transcripts)
#'     \item \code{coefficient}: Numeric coefficients associated with each feature
#'   }
#' @param plots Optional list of plots associated with the signature (default: NULL)
#' @param params Optional list of parameters used to generate the signature (default: NULL)
#' @param metadata Optional list of metadata (default: NULL)
#' @param call Optional function call (default: NULL)
#'
#' @return
#' An object of class \code{"reboot_signature"} containing:
#' \itemize{
#'   \item \code{signature}: Cleaned and sorted data.frame
#'   \item \code{plots}: Optional list of plots
#'   \item \code{params}: Optional parameter list
#'   \item \code{metadata}: Optional metadata list
#'   \item \code{call}: Function call
#' }
#'
#' @keywords internal
new_reboot_signature <- function(signature, plots = NULL, params = NULL, metadata = NULL, call = NULL) {
  
  # Checks if Reboot signature is a data frame
  if (!is.data.frame(signature)) {stop("'signature' must be a data.frame")}
  
  # Checks if signature has the 2 required columns
  if (!all(c("feature", "coefficient") %in% colnames(signature))) {
    stop("'signature' must contain 'feature' and 'coefficient' columns")
  }
  
  # Removes NA coefficients
  if (any(is.na(signature$coefficient))) {signature <- signature[!is.na(signature$coefficient), , drop = FALSE]}
  
  # Validates after NA removal
  if (nrow(signature) == 0) {stop("No valid coefficients available after removing NA values")}
  
  # Ensures numeric coefficients
  if (!is.numeric(signature$coefficient)) {stop("'coefficient' column must be numeric")}
  
  # Orders by absolute coefficient (descending)
  signature <- signature[order(abs(signature$coefficient), decreasing = TRUE), , drop = FALSE]
  
  # Resets rownames
  rownames(signature) <- NULL
  
  # Build S3 object
  obj <- structure(
    list(
      signature = signature,
      plots = plots,
      params = params,
      metadata = metadata,
      call = call
    ),
    class = "reboot_signature"
  )
  
  return(obj)
}

#' Internal constructor for REBOOT survival object
#'
#' @description Creates a standardized `reboot_survival` S3 object from REBOOT survival analysis outputs.
#' This function ensures that the generated object is structurally valid and internally consistent.
#'
#' @details This constructor is used internally by \code{rebootSurvival()} to guarantee object integrity
#' and standardization across REBOOT modules.
#'
#' It performs the following steps:
#' \itemize{
#'   \item Validates mandatory components
#'   \item Ensures core objects are properly formatted
#'   \item Standardizes optional slots
#' }
#'
#' Additional components such as plots, inputs, metadata, survival objects,
#' parameters, and function call can also be attached.
#'
#' @param signature A data.frame containing the REBOOT molecular signature with at least two columns:
#'   \itemize{
#'     \item \code{feature}: Feature names (e.g., genes or transcripts)
#'     \item \code{coefficient}: Numeric coefficients associated with each feature
#'   }
#' @param score_cutoff Numeric cutoff used to define low/high REBOOT score groups.
#' @param expression_data A data.frame containing survival variables, scores, and expression data.
#' @param clinical_data Optional data.frame containing clinical covariates (default: NULL).
#' @param survival List containing univariate and optionally multivariate survival analysis outputs.
#' @param roc_result Optional ROC analysis results (default: NULL).
#' @param ph_result Proportional hazards assumption test results.
#' @param plots Optional list of plots (default: NULL).
#' @param params Optional list of parameters used during analysis (default: NULL).
#' @param metadata Optional list of metadata information (default: NULL).
#' @param call Optional function call (default: NULL).
#'
#' @return
#' An object of class \code{"reboot_survival"} containing:
#' \itemize{
#'   \item \code{signature}: REBOOT molecular signature
#'   \item \code{score_cutoff}: Numeric score cutoff
#'   \item \code{expression_data}: Expression and survival data
#'   \item \code{clinical_data}: Clinical data (if available)
#'   \item \code{survival}: Survival analysis results
#'   \item \code{roc_result}: ROC analysis results (if available)
#'   \item \code{ph_result}: Proportional hazards results
#'   \item \code{plots}: Generated plots
#'   \item \code{params}: Analysis parameters
#'   \item \code{metadata}: Metadata information
#'   \item \code{call}: Function call
#' }
#'
#' @keywords internal
new_reboot_survival <- function(signature, score_cutoff, expression_data, survival, ph_result, clinical_data = NULL,
                                roc_result = NULL, plots = NULL, params = NULL, metadata = NULL, call = NULL) {
  
  # Checks if signature is a data frame
  if (!is.data.frame(signature)) {stop("'signature' must be a data.frame")}
  
  # Checks if expression data is a data frame
  if (!is.data.frame(expression_data)) {stop("'expression_data' must be a data.frame")}
  
  # Checks if clinical data is valid
  if (!is.null(clinical_data) && !is.data.frame(clinical_data)) {stop("'clinical_data' must be a data.frame or NULL")}
  
  # Checks if score cutoff is valid
  if (!is.numeric(score_cutoff) || length(score_cutoff) != 1) {stop("'score_cutoff' must be a single numeric value")}
  
  # Checks if survival object is valid
  if (!is.list(survival)) {stop("'survival' must be a list")}
  
  # Checks if PH result is valid
  if (is.null(ph_result)) {stop("'ph_result' cannot be NULL")}
  
  # Build S3 object
  obj <- structure(
    list(
      signature = signature,
      score_cutoff = score_cutoff,
      expression_data = expression_data,
      clinical_data = clinical_data,
      survival = survival,
      roc_result = roc_result,
      ph_result = ph_result,
      plots = plots,
      params = params,
      metadata = metadata,
      call = call
    ),
    class = "reboot_survival"
  )
  
  return(obj)
}

#' Create a REBOOT complete analysis object
#'
#' @description Internal constructor for objects of class `reboot_complete`.
#'
#' @param regression Object of class `reboot_signature` generated by \code{rebootRegression()}.
#' @param survival Object of class `reboot_survival` generated by \code{rebootSurvival()}.
#' @param params List containing parameters used in the complete pipeline.
#' @param metadata List containing metadata information.
#' @param call Original function call.
#'
#' @details Validates and creates a standardized \code{"reboot_complete"} object used internally by the package.
#'
#' @return An object of class \code{"reboot_complete"}.
#'
#' @keywords internal
new_reboot_complete <- function(regression, survival, params = NULL, metadata = NULL, call = NULL) {
  
  # Validates regression object
  if (!inherits(regression, "reboot_signature")) {stop("'regression' must be a 'reboot_signature' object")}
  
  # Validates survival object
  if (!inherits(survival, "reboot_survival")) {stop("'survival' must be a 'reboot_survival' object")}
  
  # Validates params and metadata
  if (!is.null(params) && !is.list(params)) {stop("'params' must be NULL or a list")}
  if (!is.null(metadata) && !is.list(metadata)) {stop("'metadata' must be NULL or a list")}
  
  # Creates object
  structure(
    list(
      regression = regression,
      survival = survival,
      params = params,
      metadata = metadata,
      call = call
    ),
    class = "reboot_complete"
  )
}
