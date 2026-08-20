#' Run REBOOT pipeline - module 1 (regression)
#'
#' @description
#' Main function for finding a molecular signature using bootstrap resampling and penalized (LASSO) Cox regression.
#'
#' @param filein Character or data.frame. A data.frame or a path to a tab-separated file containing survival data
#' (e.g., OS, OS.time) and expression values (genes or transcripts). Example: \code{"myPATH/myTPM.tsv"}.
#' @param outprefix Character. Output prefix used for generated files (Default: \code{"reboot"}).
#' @param bootstrap Integer. Number of bootstrap iterations (Default: \code{1}).
#' @param groupsize Integer. Number of genes/transcripts sampled per iteration (Default: \code{10}).
#' @param percentagefilter Numeric. Correlation filter threshold (0–1; Default: \code{0.3}).
#' @param variancefilter Numeric. Variance filter threshold (0–1; Default: \code{0.01}).
#' @param followup Numeric. Maximum followup time (Default: \code{NULL}).
#' @param type Character. Feature type: either \code{"gene"} (default) or \code{"transcript"}.
#' @param force Logical. If \code{TRUE}, overrides filtering steps related to survival assumptions (Default: \code{FALSE}).
#' @param plots Logical. If \code{TRUE}, generates histogram and lollipop plots (Default: \code{FALSE}).
#'
#'   These can also be generated later using \code{reboot_histogram()} and \code{reboot_lollipop()}.
#' @param table Logical. If \code{TRUE}, writes signature table to disk (Default: \code{FALSE}).
#'
#'   This can also be done later using \code{write_reboot_signature()} for the signature file.
#' @param saveJSON Logical. If \code{TRUE}, saves to disk the metadata in a JSON file (Default: \code{FALSE}).
#'
#'   This can also be done later using \code{write_reboot_metadata()} for a JSON with metadata.
#' @param saveRDS Logical. If \code{TRUE}, saves to disk the R object in a RDS file (Default: \code{FALSE}).
#'
#'   This can also be done later using \code{write_reboot_rds()}.
#' @param report Logical or character. Whether to generate a final analysis report. Accepted values are:
#' \itemize{
#'   \item \code{FALSE}: Do not generate a report (default).
#'   \item \code{TRUE}: Generate a PDF report.
#'   \item \code{"PDF"}: Generate a PDF report.
#'   \item \code{"HTML"}: Generate an HTML report.
#' }
#'
#'   This can also be done later using \code{write_reboot_report()}.
#' @param log Logical. If \code{TRUE}, writes a minimalist log file tracking the pipeline execution (Default: \code{FALSE}).
#' @param seed Integer or \code{NULL}. Random seed controlling every stochastic step in the pipeline
#'   (missing-value imputation, bootstrap feature subsampling, and LASSO cross-validation folds), so
#'   that results are reproducible across runs for the same input and parameters (Default: \code{17}).
#'
#'   Set to \code{NULL} to disable seeding and restore the previous (non-reproducible) behavior.
#'   The seed is applied locally and the caller's R session RNG state is restored afterwards.
#'
#' @details
#' The REBOOT pipeline applies a sequence of preprocessing and modeling steps:
#' \enumerate{
#'   \item Checking minimum requirements in input dataset and parameters
#'   \item Max followup definition to control for survival analysis
#'   \item Variance filtering to remove low-variance features
#'   \item Proportional hazards assumption testing (Schoenfeld residuals)
#'   \item Optional simplified regression (if dataset is small)
#'   \item Bootstrap resampling with LASSO Cox regression
#'   \item Features (pairs) correlation filtering for every iteration
#'   \item Custom coefficient filtering for genes/transcripts
#'   \item Aggregation and ranking of selected features
#' }
#'
#' The final signature consists of features with stable and non-zero coefficients across bootstrap iterations.
#'
#' @return
#' An S3 object of class \code{"reboot_signature"} containing:
#' \itemize{
#'   \item \code{signature}: data.frame with selected features and coefficients
#'   \item \code{plots}: list of generated plots (if requested)
#'   \item \code{params}: list of parameters used in the analysis
#'   \item \code{metadata}: list of metadata information such as timestamp, R version, and REBOOT package version
#'   \item \code{call}: original function call
#' }
#'
#' @note
#' Runtime may increase substantially with larger values of \code{bootstrap} (e.g., >1000 iterations) in combination with \code{groupsize}.
#'   Large dataset sizes (samples and/or features) may also significantly increase runtime.
#'
#'   For exploratory analyses, smaller values are recommended.
#'   Input data must be properly formatted and free of missing or infinite values after preprocessing steps to avoid errors.
#'
#' Change \code{percentagefilter} and \code{variancefilter} parameters to make the filters more or less flexible and,
#' as a result, obtaining bigger or smaller signatures.
#'
#' Change \code{followup} to adjust the maximum followup. Ideal when time-specific survival analyses are needed.
#' For instance, 1-year OS, or 3-year OS, or 5-year OS or even 10-year OS.
#'
#' Careful with the \code{type} parameter to get proper results, depending on the input (genes or transcripts).
#'
#' Keeping \code{force = FALSE} parameter is recommended, especially if other endpoints can be evaluated.
#'
#' @seealso
#' Project Github repository: \url{https://github.com/galantelab/reboot} \cr
#'
#' Documentation and tutorials: \url{https://galantelab.github.io/reboot/} \cr
#'
#' Other functions: \code{\link{reboot_histogram}}, \code{\link{reboot_lollipop}}, \code{\link{read_reboot_table}},
#' \code{\link{read_reboot_rds}}, \code{\link{write_reboot_rds}}, \code{\link{write_reboot_metadata}},
#' \code{\link{write_reboot_signature}}, and \code{\link{write_reboot_report}}
#'
#' @examples
#' \dontrun{
#' # Example usage
#' result <- rebootRegression(
#'   filein = "expression_data.tsv",
#'   outprefix = "my_analysis",
#'   bootstrap = 100,
#'   groupsize = 10,
#'   percentagefilter = 0.3,
#'   variancefilter = 0.01,
#'   followup = NULL,
#'   type = "gene",
#'   force = TRUE,
#'   plots = TRUE,
#'   table = TRUE,
#'   saveJSON = TRUE,
#'   saveRDS = TRUE,
#'   report = TRUE,
#'   log = TRUE
#' )
#'
#' # Access signature through custom REBOOT S3 methods
#' as.data.frame(result)
#' summary(result)
#' print(result)
#' coef(result)
#' plot(result)
#'
#' # Load previous result saved in RDS
#' read_reboot_rds(file = "my_analysis.rds")
#'
#' # Import expression and survival data manually
#' full_data <- read_reboot_table("expression_data.tsv", row.names = TRUE)
#'
#' # Generate plots later
#' reboot_histogram(result$signature, "my_analysis_hist")
#' reboot_lollipop(result$signature, "my_analysis_lolli")
#'
#' # Export signature manually
#' write_reboot_signature(result, "my_analysis_signature.txt")
#' 
#' # Export metadata manually
#' write_reboot_metadata(module = "rebootRegression", outprefix = outprefix,
#'                       parameters = list(...), command = "none")
#'
#' # Save result in RDS
#' write_reboot_rds(object = result, file = "my_analysis.rds")
#'
#' # Generate final report
#' write_reboot_report(object = result, format = report, file = "my_analysis_report.pdf")
#' }
#'
#' @export
rebootRegression <- function(filein,
                             outprefix = "reboot",
                             bootstrap = 1,
                             groupsize = 10,
                             percentagefilter = 0.3,
                             variancefilter = 0.01,
                             followup = NULL,
                             type = "gene",
                             force = FALSE,
                             plots = FALSE,
                             table = FALSE,
                             saveJSON = FALSE,
                             saveRDS = FALSE,
                             report = FALSE,
                             log = FALSE,
                             seed = 17)
{
  # Start time counter
  start_time <- Sys.time()
  
  # Checks LOG argument
  if (!is.logical(log) || length(log) != 1) {log_stop("Argument 'log' must be a boolean [TRUE or FALSE]")}
  
  # Checks SEED argument
  .reboot_check_seed(seed)

  # Checks output prefix argument
  if (!is.character(outprefix) || length(outprefix) != 1) {log_stop("Argument 'outprefix' must be a string")}
  
  # Logging setup - creating
  if (log) {
    logfile <- paste0(outprefix, "_regression.log")
    .reboot_env$log_con <- file(logfile, open = "wt")
    writeLines(paste0("[", .timestamp(), "] [INFO] REBOOT regression log started"), .reboot_env$log_con)
  }
  
  # Logging setup - closing
  on.exit({
    if (!is.null(.reboot_env$log_con) && isOpen(.reboot_env$log_con)) {
      writeLines(paste0("[", .timestamp(), "] [INFO] REBOOT regression finished"), .reboot_env$log_con)
      close(.reboot_env$log_con)
      .reboot_env$log_con <- NULL
    }
  }, add = TRUE)
  
  # Sets a local, reproducible seed covering every stochastic step used downstream: missing-value
  # imputation (mice::mice), bootstrap feature subsampling, and LASSO Cox regression's internal
  # cross-validation fold assignment (penalized::profL1). The previous global RNG state (if any)
  # is restored on exit, so this never leaks into the caller's R session. Set 'seed = NULL' to
  # fall back to ordinary (non-reproducible) random behavior.
  if (!is.null(seed)) {
    .reboot_save_seed()
    on.exit(.reboot_restore_seed(), add = TRUE)
    set.seed(seed)
  }

  # Handles 2 scenarios: (I) file PATH; (II) R data.frame
  if (is.character(filein))
  {
    if (!file.exists(filein)) {log_stop("Input file does not exist: ", filein)}
  } else if (!is.data.frame(filein)) {log_stop("'filein' must be either a valid file path or a data.frame")}
  
  # Validates input parameters
  if (!log) {log_message("REBOOT regression log started")}
  log_message("Validating parameters...")
  if (!type %in% c("gene", "transcript")) {log_stop("Argument 'type' must be either 'gene' or 'transcript'")}
  if (!is.numeric(bootstrap) || length(bootstrap) != 1 || bootstrap < 0 || bootstrap %% 1 != 0)
  {
    log_stop("Argument 'bootstrap' must be an integer")
  }
  if (!is.numeric(groupsize) || length(groupsize) != 1 || groupsize < 0 || groupsize %% 1 != 0)
  {
    log_stop("Argument 'groupsize' must be an integer")
  }
  if (!is.numeric(percentagefilter) || length(percentagefilter) != 1 ||
      percentagefilter <= 0 || percentagefilter > 1)
  {
    log_stop("Argument 'percentagefilter' must be a positive numeric value in range (0, 1]")
  }
  if (!is.numeric(variancefilter) || length(variancefilter) != 1 ||
      variancefilter <= 0 || variancefilter > 1)
  {
    log_stop("Argument 'variancefilter' must be a positive numeric value in range (0, 1]")
  }
  if (!is.null(followup) && (!is.numeric(followup) || length(followup) != 1 || followup <= 0))
  {
    log_stop("Argument 'followup' must be a positive numeric value")
  }
  if (!is.logical(force) || length(force) != 1) {log_stop("Argument 'force' must be a boolean [TRUE or FALSE]")}
  if (!is.logical(table) || length(table) != 1) {log_stop("Argument 'table' must be a boolean [TRUE or FALSE]")}
  if (!is.logical(saveJSON) || length(saveJSON) != 1) {log_stop("Argument 'saveJSON' must be a boolean [TRUE or FALSE]")}
  if (!is.logical(saveRDS) || length(saveRDS) != 1) {log_stop("Argument 'saveRDS' must be a boolean [TRUE or FALSE]")}
  if (!is.logical(plots) || length(plots) != 1) {log_stop("Argument 'plots' must be a boolean [TRUE or FALSE]")}
  if (!(is.logical(report) || (is.character(report) && length(report) == 1 && toupper(report) %in% c("PDF", "HTML"))))
  {
    log_stop("'report' must be a boolean [TRUE or FALSE] or a character ['PDF' or 'HTML']")
  }
  log_message("Done.")
  
  # Reads INput
  log_message("Reading input expression/survival file...")
  if (is.character(filein)) {full_data <- read_reboot_table(filein, sep = "\t")} else {full_data <- filein}
  colnames(full_data) <- gsub("-", "__", colnames(full_data))
  log_message("Done.")
  
  # Standardizes the column names of the survival metrics
  log_message("Renaming first two (survival) columns to 'OS' and 'OS.time'...")
  colnames(full_data)[1:2] <- c("OS", "OS.time")
  log_message("Done.")
  
  # Validates input file
  log_message("Validating input...")
  
  # General structure
  if (ncol(full_data) < 3) {log_stop("Input data must contain at least two [OS, OS.time] survival columns and one feature")}
  
  # OS survival column
  os_values <- full_data[,1]
  if (!is.numeric(os_values)) {log_stop("Column [OS] must be numeric and the first in the table")}
  if (any(is.na(os_values))) {log_stop("Column [OS] must be the first in the table and contain no missing values")}
  unique_os <- sort(unique(stats::na.omit(os_values)))
  if (!setequal(unique_os, c(0, 1))) {log_stop("Column [OS] must be the first in the table and contain only binary values: 0 and 1")}
  
  # OS.time survival column
  time_values <- full_data[,2]
  if (!is.numeric(time_values)) {log_stop("Column [OS.time] must be numeric and the second in the table")}
  if (any(is.na(time_values))) {log_stop("Column [OS.time] must be the second in the table and contain no missing values")}
  if (any(time_values < 0, na.rm = TRUE)) {log_stop("Column [OS.time] must be the second in the table and contain only positive values")}
  log_message("Done.")
  
  # Applies maximum follow-up censoring if requested
  log_message("Applying maximum followup value...")
  full_data <- reboot_followup(data = full_data, followup = followup)
  log_message("Done.")
  
  # Defines coefficient cutoffs to minimize false positives
  bar <- if (type == "gene") 0.0035 else 0.011
  
  # Variance filter
  log_message("Running variance filter...")
  full_data <- reboot_varfun(full_data, variancefilter, force)
  log_message("Done.")
  
  # Schoenfeld tests
  log_message("Performing Schoenfeld test...")
  full_data <- reboot_ph_assumptions(full_data)
  log_message("Done.")
  
  # Validates processed data
  log_message("Validating processed data...")
  if (!is.data.frame(full_data)) {log_stop("Processed data is not a 'data.frame' anymore")}
  if (nrow(full_data) == 0) {log_stop("No samples remaining after filtering steps")}
  if ((ncol(full_data) - 2) < 1) {log_stop("No features remaining after filtering steps")}
  num_data <- full_data[, sapply(full_data, is.numeric), drop = FALSE]
  if (any(!is.finite(as.matrix(num_data)))) {log_stop("Processed data contains invalid values [NA, NaN or Inf]")}
  log_message("Done.")
  
  # Runs simplified regression
  log_message("Checking if simplified regression with no bootstraps is required...")
  early_signature <- reboot_feature_check(full_data, groupsize)
  log_message("Done.")
  
  # Decides if user needs to use simplified result or go for the more robust Bootstrap approach
  if (!is.null(early_signature)) {
    log_message("Using simplified regression with no bootstraps...")
    signature <- early_signature
    log_message("Done.")
  } else {
    log_message("Running robust regression with bootstraps...")
    signature <- reboot_bootstrapfun(full_data, bootstrap, groupsize, percentagefilter, bar)
    log_message("Done.")
  }
  
  # Validation
  log_message("Validating generated signature...")
  if (!is.data.frame(signature)) {log_stop("signature must be a 'data.frame'")}
  if (nrow(signature) == 0) {log_stop("No features selected in the final signature")}
  if (!all(c("feature", "coefficient") %in% colnames(signature))) {log_stop("signature must contain 'feature' and 'coefficient' columns")}
  if (any(is.na(signature$coefficient)))
  {
    log_warning("NA values found in coefficients. Removing them...")
    signature <- signature[!is.na(signature$coefficient), ]
  }
  if (nrow(signature) == 0) {log_stop("All coefficients were NA after filtering")}
  signature <- signature[order(abs(signature$coefficient), decreasing = TRUE), ]
  log_message("Done.")
  
  # Produces histogram and lollipop plots for generated signature
  plot_list <- NULL
  if (plots)
  {
    log_message("Making histogram and lollipop plots...")
    plot_list <- list(
      hist_plot = reboot_histogram(signature, outprefix),
      lolli_plot = reboot_lollipop(signature, outprefix))
    log_message("Done.")
  }
  
  # Writes signature table if requested
  if (table)
  {
    log_message("Writing signature table...")
    write_reboot_signature(signature, paste0(outprefix, "_signature.txt"))
    log_message("Done.")
  }
  
  # Creates S3 object containing the signature and plots
  # Listed parameters are only those directly affecting the scientific results/model
  result <- new_reboot_signature(
    signature = signature,
    plots = plot_list,
    params = list(
      bootstrap = bootstrap,
      groupsize = groupsize,
      percentagefilter = percentagefilter,
      variancefilter = variancefilter,
      followup = followup,
      type = type,
      force = force,
      seed = seed
    ),
    metadata = list(
      timestamp = Sys.time(),
      R_version = R.version.string,
      package_version = utils::packageVersion("Reboot")
    ),
    call = match.call()
  )
  
  # Writes execution metadata
  if (saveJSON)
  {
    log_message("Writing JSON metadata...")
    write_reboot_metadata(
      module = "rebootRegression",
      outprefix = outprefix,
      parameters = list(
        filein = if (is.character(filein)) {filein} else {"data.frame"},
        outprefix = outprefix,
        bootstrap = bootstrap,
        groupsize = groupsize,
        percentagefilter = percentagefilter,
        variancefilter = variancefilter,
        followup = followup,
        type = type,
        force = force,
        plots = plots,
        table = table,
        saveJSON = saveJSON,
        saveRDS = saveRDS,
        report = report,
        log = log,
        seed = seed
      ),
      command = "none"
    )
    log_message("Done.")
  }
  
  # Writes RDS object
  if (saveRDS)
  {
    log_message("Writing RDS object...")
    write_reboot_rds(object = result, file = paste0(outprefix, "_regression.rds"))
    log_message("Done.")
  }
  
  # Generates report if requested
  if (!isFALSE(report))
  {
    log_message("Writing REBOOT report...")
    write_reboot_report(object = result, format = report, file = paste0(outprefix, "_report_regression.pdf"))
    log_message("Done.")
  }
  
  # End time counter
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Prints elapsed time
  if (elapsed_time >= 3600) {
    log_message(sprintf("Time to run [module I]: %.2f hours", elapsed_time / 3600))
  } else if (elapsed_time >= 60) {
    log_message(sprintf("Time to run [module I]: %.2f minutes", elapsed_time / 60))
  } else {
    log_message(sprintf("Time to run [module I]: %.2f seconds", elapsed_time))
  }
  
  # Returns the OBJ with the signature
  if (!log) {log_message("REBOOT regression finished")}
  return(result)
}

#' Run REBOOT pipeline - module 2 (survival analysis)
#'
#' @description
#' Applies a molecular signature to survival data, generating prognostic scores and performing
#' survival analyses (univariate and, optionally, multivariate Cox regression).
#'
#' @param filein Character or data.frame. A data.frame or a path to a tab-separated file containing survival data
#' (e.g., OS, OS.time) and expression values (genes or transcripts). Example: \code{"myPATH/myTPM.tsv"}.
#' @param signature A molecular signature. Can be:
#'   \itemize{
#'     \item An object of class \code{"reboot_signature"} from the \code{rebootRegression()} Reboot function (regression - module I)
#'     \item A data.frame with columns \code{feature} and \code{coefficient}
#'     \item A path to a .tsv file containing the signature
#'   }
#' @param outprefix Character. Output prefix used for generated files (Default: \code{"reboot"}).
#' @param multivariate Logical. If \code{TRUE}, performs multivariate survival analysis including clinical variables
#'   present in \code{clinin} (Default: \code{FALSE}).
#' @param clinin Character or data.frame. A data.frame or a path to a tab-separated file containing binary-categorized
#' clinical data (e.g., Age, Gender, ...) (required if \code{multivariate = TRUE}). Example: \code{"myPATH/myClinics.tsv"}.
#' @param roc Logical. If \code{TRUE}, uses a ROC curve to define the score cutoff instead of median (Default: \code{FALSE}).
#' @param variancefilter Numeric. Minimum normalized variance for follow-up time (0–1; Default: \code{0.01}).
#' @param followup Numeric. Maximum followup time (Default: \code{NULL}).
#' @param p.cutoff Numeric. P-value threshold used to select clinical covariates from univariate Cox regression
#' for inclusion in multivariate analysis (Default: \code{0.2}).
#' @param bootstrap Integer. Number of bootstrap iterations for multivariate Cox regression (Default: \code{1}).
#' @param force Logical. If \code{TRUE}, overrides filtering steps related to survival assumptions (Default: \code{FALSE}).
#' @param plots Logical. If \code{TRUE}, generates a Kaplan-Meier curve, the proportional hazards assumption plot, and
#'   extra figures, depending on the chosen parameters (Default: \code{FALSE}).
#' 
#'   These can also be generated later using \code{reboot_km()} and \code{reboot_phplot()}.
#'   If \code{roc = TRUE}, a ROC curve can also be generated using \code{reboot_roc()}.
#'   If \code{multivariate = TRUE}, a forest plot can also be generated using \code{reboot_forest()}.
#'   If bootstrap is performed, an extra summary barplot can also be generated using \code{reboot_barplot()}.
#' @param table Logical. If \code{TRUE}, writes output tables to disk (Default: \code{FALSE}).
#'
#'   These can also be done later using \code{write_reboot_signature()} for an updated signature file,
#'   \code{write_reboot_score()} for a table showing the calculated Reboot scores, \code{write_reboot_univariate()}
#'   for a file containing the univariate survival results, and \code{write_reboot_multivariate()}
#'   for a file containing the multivariate survival results.
#' @param saveJSON Logical. If \code{TRUE}, saves to disk the metadata in a JSON file (Default: \code{FALSE}).
#'
#'   This can also be done later using \code{write_reboot_metadata()} for a JSON with metadata.
#' @param saveRDS Logical. If \code{TRUE}, saves to disk the R object in a RDS file (Default: \code{FALSE}).
#'
#'   This can also be done later using \code{write_reboot_rds()} for a RDS object.
#' @param report Logical or character. Whether to generate a final analysis report. Accepted values are:
#' \itemize{
#'   \item \code{FALSE}: Do not generate a report (default).
#'   \item \code{TRUE}: Generate a PDF report.
#'   \item \code{"PDF"}: Generate a PDF report.
#'   \item \code{"HTML"}: Generate an HTML report.
#' }
#'
#'   This can also be done later using \code{write_reboot_report()}.
#' @param log Logical. If \code{TRUE}, writes a minimalist log file tracking the pipeline execution (Default: \code{FALSE}).
#' @param seed Integer or \code{NULL}. Random seed controlling the stochastic bootstrap resampling
#'   used to select stable clinical covariates in multivariate analysis, so that results are
#'   reproducible across runs for the same input and parameters (Default: \code{17}).
#'
#'   Set to \code{NULL} to disable seeding and restore the previous (non-reproducible) behavior.
#'   The seed is applied locally and the caller's R session RNG state is restored afterwards.
#'
#' @details
#' The REBOOT pipeline applies a sequence of analyses steps:
#' \enumerate{
#'   \item Checking minimum requirements in input datasets and parameters
#'   \item Max followup definition to control for survival analysis
#'   \item Updating the molecular signature
#'   \item Creating the signature scores for each patient/sample
#'   \item Proportional hazards assumption testing (Schoenfeld residuals)
#'   \item Categorizing the scores into low/high groups
#'   \item Performing the univariate Cox regression for the score
#'   \item Running multiple univariate Cox regression for other clinical variables
#'   \item Performing the multivariate Cox regression with bootstrap resampling if the requirements are met
#' }
#'
#' The final complete results consist of multiple tables and figures to show the association between the Reboot scores,
#' created based on the molecular signature, and survival parameters, optionally correcting for clinical variables.
#'
#' @return
#' An S3 object of class \code{"reboot_survival"} containing:
#' \itemize{
#'   \item \code{signature}: data.frame with selected features, coefficients, and extra information
#'   \item \code{score_cutoff}: cutoff used to create the low/high score groups
#'   \item \code{expression_data}: data.frame with scores and survival information for all patients/samples
#'   \item \code{clinical_data}: data.frame with score groups and clinical variables for all patients/samples (if requested)
#'   \item \code{survival}: survival analysis outputs
#'   \item \code{roc_result}: ROC analysis outputs (if requested)
#'   \item \code{ph_result}: proportional hazards analysis outputs
#'   \item \code{plots}: list of generated plots (if requested)
#'   \item \code{params}: list of parameters used in the analysis
#'   \item \code{metadata}: list of metadata information such as timestamp, R version, and REBOOT package version
#'   \item \code{call}: original function call
#' }
#'
#' @note
#' Runtime may increase substantially with larger values of \code{bootstrap} (e.g., >1000 iterations) and
#' large dataset sizes (samples and/or features and clinical variables).
#'
#'   For exploratory analyses, smaller values are recommended.
#'   Input data must be properly formatted and free of missing or infinite values after preprocessing steps to avoid errors.
#'
#' Users are allowed and encouraged to further investigate the Reboot \code{signature} provided. One possibility is to
#' test smaller signatures (e.g., establishing a coefficient cutoff) and even to find a minimum signature related to
#' patient survival.
#'
#' Change \code{followup} to adjust the maximum followup. Ideal when time-specific survival analyses are needed.
#' For instance, 1-year OS, or 3-year OS, or 5-year OS or even 10-year OS.
#'
#' Keeping \code{force = FALSE} parameter is recommended, especially if other endpoints can be evaluated.
#' In case an error occurs, it is recommended to experiment changing \code{variancefilter} and \code{p.cutoff}
#' parameters first to make the follow up variance and variable selection filters more flexible.
#'
#' @seealso
#' Project Github repository: \url{https://github.com/galantelab/reboot} \cr
#'
#' Documentation and tutorials: \url{https://galantelab.github.io/reboot/} \cr
#'
#' Other functions: \code{\link{reboot_km}}, \code{\link{reboot_phplot}}, \code{\link{reboot_roc}},
#' \code{\link{reboot_forest}}, \code{\link{reboot_barplot}}, \code{\link{read_reboot_table}},
#' \code{\link{read_reboot_rds}}, \code{\link{write_reboot_rds}}, \code{\link{write_reboot_signature}},
#' \code{\link{write_reboot_score}}, \code{write_reboot_metadata()}, \code{\link{write_reboot_univariate}},
#' \code{\link{write_reboot_multivariate}}, and \code{\link{write_reboot_report}}
#'
#' @examples
#' \dontrun{
#' # Example usage
#' result <- rebootSurvival(
#'   filein = "expression_data.tsv",
#'   signature = "signature_data.tsv",
#'   outprefix = "my_analysis",
#'   multivariate = TRUE,
#'   clinin = "clinical_data.tsv",
#'   roc = TRUE,
#'   variancefilter = 0.01,
#'   followup = NULL,
#'   p.cutoff = 0.2,
#'   bootstrap = 100,
#'   force = TRUE,
#'   plots = TRUE,
#'   table = TRUE,
#'   saveJSON = TRUE,
#'   saveRDS = TRUE,
#'   report = TRUE,
#'   log = TRUE
#' )
#'
#' # Access signature through custom REBOOT S3 methods
#' as.data.frame(result)
#' summary(result)
#' print(result)
#' coef(result)
#' plot(result)
#'
#' # Load previous result saved in RDS
#' read_reboot_rds(file = "my_analysis.rds")
#'
#' # Import signature data manually
#' raw_sig_df <- read_reboot_table("signature_data.tsv", row.names = FALSE)
#'
#' # Import expression and survival data manually
#' data <- read_reboot_table("expression_data.tsv", row.names = TRUE)
#'
#' # Import clinical data manually
#' clin <- read_reboot_table("clinical_data.tsv", row.names = TRUE)
#'
#' # Generate plots later (assuming multivariate analysis was performed)
#' reboot_roc(x = result$roc_result, outprefix = "my_analysis_roc")
#' reboot_km(data = result$expression_data, cutoff = result$score_cutoff,
#'           pval = result$survival$univariate$table$log.rank.pvalue, outprefix = "my_analysis_km")
#' reboot_phplot(x = result$ph_result, is_multi = TRUE, outprefix = "my_analysis_phplot")
#' reboot_forest(model_object = result$survival$multivariate$model, outprefix = "my_analysis_forest")
#' reboot_barplot(plot_df = result$survival$multivariate$bootstrap.table,
#'                outprefix = "my_analysis_bar")
#'
#' # Export signature manually
#' write_reboot_signature(result, "my_analysis_signature_updated.tsv")
#'
#' # Export scores manually
#' write_reboot_score(result, "my_analysis_scoreCont.tsv") # for univariate mode only
#' write_reboot_score(result, "my_analysis_scoreCat.tsv", clinics = TRUE) # if multivariate mode is on
#'
#' # Export survival results manually
#' write_reboot_univariate(result, "my_analysis_uniCox.txt") # for univariate mode only
#' write_reboot_multivariate(result, "my_analysis_multiCox.txt") # if multivariate mode is on
#'
#' # Export metadata manually
#' write_reboot_metadata(module = "rebootSurvival", outprefix = outprefix,
#'                       parameters = list(...), command = "none")
#'
#' # Save result in RDS
#' write_reboot_rds(object = result, file = "my_analysis.rds")
#'
#' # Generate final report
#' write_reboot_report(object = result, format = report, file = "my_analysis_report.pdf")
#' }
#'
#' @export
rebootSurvival <- function(filein,
                           signature,
                           outprefix = "reboot",
                           multivariate = FALSE,
                           clinin = NULL,
                           roc = FALSE,
                           variancefilter = 0.01,
                           followup = NULL,
                           p.cutoff = 0.2,
                           bootstrap = 1,
                           force = FALSE,
                           plots = FALSE,
                           table = FALSE,
                           saveJSON = FALSE,
                           saveRDS = FALSE,
                           report = FALSE,
                           log = FALSE,
                           seed = 17)
{
  # Start time counter
  start_time <- Sys.time()
  
  # Initializes specific multivariate objects
  multi_cox <- NULL
  multi_model <- NULL
  clin <- NULL
  
  # Checks LOG argument
  if (!is.logical(log) || length(log) != 1) {log_stop("Argument 'log' must be a boolean [TRUE or FALSE]")}
  
  # Checks SEED argument
  .reboot_check_seed(seed)

  # Checks output prefix argument
  if (!is.character(outprefix) || length(outprefix) != 1) {log_stop("Argument 'outprefix' must be a string")}
  
  # Logging setup - creating
  if (log) {
    logfile <- paste0(outprefix, "_survival.log")
    .reboot_env$log_con <- file(logfile, open = "wt")
    writeLines(paste0("[", .timestamp(), "] [INFO] REBOOT survival log started"), .reboot_env$log_con)
  }
  
  # Logging setup - closing
  on.exit({
    if (!is.null(.reboot_env$log_con) && isOpen(.reboot_env$log_con)) {
      writeLines(paste0("[", .timestamp(), "] [INFO] REBOOT survival finished"), .reboot_env$log_con)
      close(.reboot_env$log_con)
      .reboot_env$log_con <- NULL
    }
  }, add = TRUE)
  
  # Sets a local, reproducible seed covering the stochastic step used downstream: bootstrap
  # resampling of clinical covariates for the multivariate Cox model (sjstats::bootstrap, only
  # relevant when multivariate = TRUE). The previous global RNG state (if any) is restored on
  # exit, so this never leaks into the caller's R session. Set 'seed = NULL' to fall back to
  # ordinary (non-reproducible) random behavior.
  if (!is.null(seed)) {
    .reboot_save_seed()
    on.exit(.reboot_restore_seed(), add = TRUE)
    set.seed(seed)
  }

  # Handles 2 scenarios: (I) file PATH; (II) R data.frame
  if (!log) {log_message("REBOOT survival log started")}
  log_message("Validating file parameters...")
  if (is.character(filein))
  {
    if (!file.exists(filein)) {log_stop("Input file does not exist: ", filein)}
  } else if (!is.data.frame(filein)) {log_stop("'filein' must be either a valid file path or a data.frame")}
  
  # Validates file input parameters
  if (!is.logical(multivariate) || length(multivariate) != 1) {log_stop("Argument 'multivariate' must be a boolean [TRUE or FALSE]")}
  
  # Validates clinical input for multivariate analysis
  if (multivariate)
  {
    if (is.null(clinin)) {log_stop("'clinin' must be provided when 'multivariate = TRUE'")}
    if (is.character(clinin))
    {
      if (!file.exists(clinin)) {log_stop("Clinical input file does not exist: ", clinin)}
    } else if (!is.data.frame(clinin)) {log_stop("'clinin' must be either a valid file path or a data.frame")}
  }
  log_message("Done.")
  
  # Loads molecular signature
  log_message("Loading signature...")
  
  # Handles 3 scenarios: (I) signature from module I (regression); (II) user PATH; (III) R data.frame
  if (inherits(signature, "reboot_signature")) {
    raw_sig_df <- signature$signature
  } else if (is.character(signature)) {
    if (!file.exists(signature)) {log_stop("Signature file does not exist: ", signature)}
    raw_sig_df <- read_reboot_table(signature, sep = "\t", row.names = FALSE)
  } else if (is.data.frame(signature)) {
    raw_sig_df <- signature
  } else {
    log_stop("Invalid 'signature' input")
  }
  
  # Guarantees the signature is of class "reboot_signature"
  raw_sig_df <- new_reboot_signature(raw_sig_df)$signature
  log_message("Done.")
  
  # Adds extra column to associate the coefficient signal to its prognostic value
  raw_sig_df$prognostic <- ifelse(raw_sig_df$coefficient > 0, "worse", "better")
  
  # Writes updated signature table if requested
  if (table)
  {
    log_message("Writing updated signature table...")
    write_reboot_signature(raw_sig_df, paste0(outprefix, "_signature_updated.tsv"))
    log_message("Done.")
  }
  
  # Gets only the first 2 relevant columns ('feature' and 'coefficient') for downstream analyses
  sig_df <- raw_sig_df[,c(1, 2)]
  
  # Adjusts feature names to match the format of those found in the survival/expression table
  sig_df$feature <- gsub("-", "__", sig_df$feature)
  
  # Validates other parameters
  log_message("Validating other parameters...")
  if (!is.numeric(variancefilter) || length(variancefilter) != 1 || variancefilter <= 0 || variancefilter > 1)
  {
    log_stop("Argument 'variancefilter' must be a positive numeric value in range (0, 1]")
  }
  if (!is.null(followup) && (!is.numeric(followup) || length(followup) != 1 || followup <= 0))
  {
    log_stop("Argument 'followup' must be a positive numeric value")
  }
  if (!is.numeric(p.cutoff) || length(p.cutoff) != 1 || p.cutoff <= 0 || p.cutoff >= 1)
  {
    stop("'p.cutoff' must be a numeric value between 0 and 1")
  }
  if (!is.numeric(bootstrap) || length(bootstrap) != 1 || bootstrap < 0 || bootstrap %% 1 != 0)
  {
    log_stop("Argument 'bootstrap' must be an integer")
  }
  if (!is.logical(roc) || length(roc) != 1) {log_stop("Argument 'roc' must be a boolean [TRUE or FALSE]")}
  if (!is.logical(force) || length(force) != 1) {log_stop("Argument 'force' must be a boolean [TRUE or FALSE]")}
  if (!is.logical(table) || length(table) != 1) {log_stop("Argument 'table' must be a boolean [TRUE or FALSE]")}
  if (!is.logical(saveJSON) || length(saveJSON) != 1) {log_stop("Argument 'saveJSON' must be a boolean [TRUE or FALSE]")}
  if (!is.logical(saveRDS) || length(saveRDS) != 1) {log_stop("Argument 'saveRDS' must be a boolean [TRUE or FALSE]")}
  if (!is.logical(plots) || length(plots) != 1) {log_stop("Argument 'plots' must be a boolean [TRUE or FALSE]")}
  if (!(is.logical(report) || (is.character(report) && length(report) == 1 && toupper(report) %in% c("PDF", "HTML"))))
  {
    log_stop("'report' must be a boolean [TRUE or FALSE] or a character ['PDF' or 'HTML']")
  }
  log_message("Done.")
  
  # Loads INput with expression and survival data
  log_message("Reading input expression/survival file...")
  if (is.character(filein)) {data <- read_reboot_table(filein, sep = "\t")} else {data <- filein}
  colnames(data) <- gsub("-", "__", colnames(data))
  log_message("Done.")
  
  # Standardizes the column names of the survival metrics
  log_message("Renaming first two (survival) columns to 'OS' and 'OS.time'...")
  colnames(data)[1:2] <- c("OS", "OS.time")
  log_message("Done.")
  
  # Validates input file
  log_message("Validating input...")
  
  # General structure
  if (ncol(data) < 3) {log_stop("Input data must contain at least two [OS, OS.time] survival columns and one feature")}
  
  # OS survival column
  os_values <- data[,1]
  if (!is.numeric(os_values)) {log_stop("Column [OS] must be numeric and the first in the table")}
  if (any(is.na(os_values))) {log_stop("Column [OS] must be the first in the table and contain no missing values")}
  unique_os <- sort(unique(stats::na.omit(os_values)))
  if (!setequal(unique_os, c(0, 1))) {log_stop("Column [OS] must be the first in the table and contain only binary values: 0 and 1")}
  
  # OS.time survival column
  time_values <- data[,2]
  if (!is.numeric(time_values)) {log_stop("Column [OS.time] must be numeric and the second in the table")}
  if (any(is.na(time_values))) {log_stop("Column [OS.time] must be the second in the table and contain no missing values")}
  if (any(time_values < 0, na.rm = TRUE)) {log_stop("Column [OS.time] must be the second in the table and contain only positive values")}
  log_message("Done.")
  
  # Applies maximum follow-up censoring if requested
  log_message("Applying maximum followup value...")
  data <- reboot_followup(data = data, followup = followup)
  log_message("Done.")
  
  # Applies OS and OS.time filters if user set force = FALSE (default)
  if (!force)
  {
    perc <- sum(data$OS) / nrow(data)
    if (perc < 0.2 || perc > 0.8)
    {
      log_stop(paste0("Survival status proportion outside acceptable range [between 20-80%]. ",
                      "Consider other endpoints. If you want to continue anyway, set force = T."))
    }
    fvar <- stats::var(data$OS.time / max(data$OS.time))
    if (fvar < variancefilter)
    {
      log_stop(paste0("Follow-up variance: ", fvar, " below threshold. Consider changing the 'variancefilter' parameter",
                      " or choosing other endpoints. If you want to continue anyway, set force = T."))
    }
  }
  
  # Checks if all features in signature are present in the expression file
  log_message("Checking compatibility between expression/survival and signature files...")
  expr_features <- colnames(data)[-(1:2)]
  if (anyDuplicated(expr_features)) {log_stop("Duplicated features are not allowed in expression data")}
  if (anyDuplicated(sig_df$feature)) {log_stop("Duplicated features are not allowed in signature")}
  missing_features <- setdiff(sig_df$feature, expr_features)
  if (length(missing_features) > 0) {log_stop("Signature features missing in expression data: ", paste(missing_features, collapse = ", "))}
  log_message("Done.")
  
  # Reorders and subsets data according to signature features and order, essential to correct score calculation later
  data <- data[, c("OS", "OS.time", sig_df$feature), drop = FALSE]
  
  # Calculates Reboot score
  log_message("Calculating signature score...")
  
  # Reorders expression matrix according to signature features
  expr <- as.matrix(data[, sig_df$feature, drop = FALSE])
  
  # Explicitly aligns coefficients to expression columns
  match_idx <- match(colnames(expr), sig_df$feature)
  
  # Safety check
  if (any(is.na(match_idx))) {log_stop("Expression matrix and signature coefficients are not aligned")}
  coef_vector <- sig_df$coefficient[match_idx]
  
  # Calculates signature score
  score <- expr %*% coef_vector
  
  # Ensures the calculated scores are numeric
  data$score <- as.numeric(score)
  log_message("Done.")
  
  # Categorizes the Reboot score into high/low groups
  log_message("Categorizing score...")
  
  # Based on user's choice, use the ROC curve or the median
  if (roc)
  {
    # This value can be changed. By default, it uses the median followup time
    roc_obj <- survivalROC::survivalROC(Stime = data$OS.time, status = data$OS, marker = data$score,
                                        predict.time = stats::median(data$OS.time, na.rm = TRUE), method = "NNE",
                                        span = 0.25*nrow(data)^(-0.20))
    roc_cutoff <- reboot_cutoff_roc(data, direction = ifelse(as.numeric(roc_obj$AUC) < 0.5, ">", "<"))
    score_cutoff <- roc_cutoff$cutoff
    log_message(paste0("ROC cutoff: ", score_cutoff))
  } else {
    score_cutoff <- stats::median(data$score, na.rm = TRUE)
    log_message(paste0("Median cutoff: ", score_cutoff))
  }
  
  # Explicitly set the "high" group as REF to control for HR (CI 95%) values
  data$score_group <- ifelse(data$score > score_cutoff, "high", "low")
  data$score_group <- factor(data$score_group, levels = c("high", "low"))
  log_message("Done.")
  
  # Writes table with scores (both values and categories) if requested
  if (table)
  {
    log_message("Writing score table...")
    colnames(data)[colnames(data) == "score_group"] <- "Score"
    write_reboot_score(x = data, file = paste0(outprefix, "_scoreCont.tsv"), clinics = FALSE)
    colnames(data)[colnames(data) == "Score"] <- "score_group"
    log_message("Done.")
  }
  
  # Performs univariate survival analysis
  log_message("Running univariate survival analysis...")
  
  # Runs log-rank test for score signatures
  res_logrank <- reboot_uniCox_model(data)
  log_message("Done.")
  
  # Checks if model is valid
  if (is.null(res_logrank$table)) {log_stop("Cox model failed (convergence issue)")}
  
  # Writes table with univariate survival results if requested
  if (table)
  {
    log_message("Writing univariate survival results table...")
    res_logrank$table$feature[res_logrank$table$feature == "score"] <- "Score"
    write_reboot_univariate(res_logrank$table, paste0(outprefix, "_uniCox.txt"))
    res_logrank$table$feature[res_logrank$table$feature == "Score"] <- "score"
    log_message("Done.")
  }
  
  # Performs multivariate survival analysis
  if (multivariate)
  {
    log_message("Running multivariate survival analysis...")
    
    # Reads clinical data
    log_message("Reading input clinical file...")
    if (is.character(clinin)) {clin <- read_reboot_table(clinin, sep = "\t")} else {clin <- clinin}
    log_message("Done.")
    
    # Checks if all samples are present in both datasets
    log_message("Checking compatibility between expression/survival and clinical files...")
    if (!identical(sort(rownames(clin)), sort(rownames(data))))
    {
      miss_in_clin <- setdiff(rownames(data), rownames(clin))
      miss_in_data <- setdiff(rownames(clin), rownames(data))
      if (length(miss_in_clin) > 0) {log_stop("Samples missing in clinical data: ", paste(miss_in_clin, collapse = ", "))}
      if (length(miss_in_data) > 0) {log_stop("Samples missing in expression/survival data: ", paste(miss_in_data, collapse = ", "))}
    }
    log_message("Done.")
    
    # Reorders clinical data according to expression/survival data
    clin <- clin[rownames(data), , drop = FALSE]
    
    # Keeps only survival and score-related columns
    data_surv <- data[, c("OS", "OS.time", "score", "score_group"), drop = FALSE]
    
    # Merges survival/score data with clinical covariates
    clin <- cbind(data_surv, clin)
    
    # Validates if all clinical variables are binary categorical
    log_message("Checking if all variables in clinical file are binary...")
    for (i in seq_len(ncol(clin))[-(1:4)])
    {
      if (nlevels(as.factor(clin[, i])) != 2)
      {
        log_stop("Clinical variable '", colnames(clin)[i], "' must contain exactly 2 categories.")
      }
    }
    log_message("Done.")
    
    # Defines clinical covariates
    log_message("Selecting covariables with multiple univariate Cox analyses...")
    uni_covariates <- colnames(clin)[5:ncol(clin)]
    
    # Runs multiple univariate Cox analyses
    univ_cox <- reboot_uniCox_test(data = clin, covariates = uni_covariates, p.cutoff = p.cutoff)
    
    # Selects covariates with p-value < user-defined cutoff (ex: 0.2)
    selected <- univ_cox[!is.na(univ_cox$Cox.pvalue) & univ_cox$Cox.pvalue < p.cutoff, , drop = FALSE]
    log_message("Done.")
    
    # Checks if any covariate passed univariate filtering
    if (nrow(selected) == 0)
    {
      log_warning("No covariables passed univariate filtering. Multivariate analysis could not be performed.")
    } else {
      # Transforms clinical covariate names
      selected_covariates <- sapply(selected$variable, function(var)
      {
        matched <- uni_covariates[sapply(uni_covariates, function(x) grepl(paste0("^", x), var))]
        if (length(matched) > 0) matched[[1]] else var
      })
      selected_covariates <- unique(selected_covariates)
      log_message("Covariables selected: ", paste(selected_covariates, collapse = ", "))
      
      # Runs initial multivariate Cox regression and builds initial model
      log_message("Building initial model...")
      multi_cox_initial <- reboot_multiCox_test(data = clin, univ_result = univ_cox, logrank_result = res_logrank$table,
                                                covariates = selected_covariates, all_covariates = uni_covariates,
                                                use.bootstrap = TRUE, bootstrap = bootstrap, bootstrap.freq = 0.25, multi.p.cutoff = 0.05)
      if (multi_cox_initial$bootstrap.used)
      {
        log_message(paste0("Bootstrap resampling performed with ", bootstrap, " iterations."))
      } else {
        if (bootstrap == 1)
        {
          log_message(paste0("Bootstrap not performed. Default value: ", bootstrap, " used."))
        } else {
          log_message("Bootstrap failed. Resuming analysis with no resampling...")
        }
      }
      log_message("Done.")
      
      # Tests proportional hazards assumptions
      log_message("Testing proportional hazards assumptions...")
      ph_test <- reboot_schoenfeld(model_object = multi_cox_initial$model, covariates = multi_cox_initial$covariates, is_multi = TRUE)
      ph_covariates <- ph_test$variables
      log_message("Done.")
      
      # Checks if any covariate passed proportional hazards assumptions
      if (length(ph_covariates) == 0 || (length(ph_covariates) == 1 && ph_covariates == "score_group"))
      {
        log_warning("No covariables passed PH assumptions. Multivariate analysis could not be performed.")
      } else {
        # Runs final multivariate Cox regression and builds final model
        log_message("Building final model...")
        multi_cox <- reboot_multiCox_test(data = clin, univ_result = univ_cox, logrank_result = res_logrank$table,
                                          covariates = ph_covariates, all_covariates = uni_covariates, use.bootstrap = FALSE)
        multi_model <- multi_cox$model
        log_message("Done.")
        
        # Writes tables with multivariate survival results and score + clinics variables if requested
        if (table)
        {
          log_message("Writing clinical and multivariate survival results tables...")
          colnames(clin)[colnames(clin) == "score_group"] <- "Score"
          write_reboot_score(x = clin, file = paste0(outprefix, "_scoreCat.tsv"), clinics = TRUE)
          colnames(clin)[colnames(clin) == "Score"] <- "score_group"
          multi_cox$table$variable[multi_cox$table$variable == "score_group"] <- "Score"
          write_reboot_multivariate(multi_cox$table, paste0(outprefix, "_multiCox.txt"))
          log_message("Done.")
        }
      }
    }
    log_message("Done.")
  }
  
  # Produces uni and multivariate plots
  plot_list <- NULL
  if (plots)
  {
    log_message("Making plots from univariate and multivariate (if applicable) survival analyses...")
    plot_list <- list(
      km_plot = reboot_km(data = data, cutoff = score_cutoff, outprefix = outprefix,
                          pval = res_logrank$table$log.rank.pvalue, palette = c("#D73027", "#1A9850")),
      ph_plot = if (!is.null(multi_model)) reboot_phplot(x = ph_test, is_multi = TRUE,
                                                         outprefix = outprefix) else reboot_phplot(x = res_logrank$ph,
                                                                                                   outprefix = outprefix),
      roc_plot = if (isTRUE(roc)) reboot_roc(x = roc_cutoff, outprefix = outprefix) else NULL,
      forest_plot = if (!is.null(multi_model)) reboot_forest(model_object = multi_model,
                                                             outprefix = outprefix) else NULL,
      bar_plot = if (!is.null(multi_cox_initial$model) &&
                     isTRUE(multi_cox_initial$bootstrap.used)) reboot_barplot(plot_df = multi_cox$bootstrap.table,
                                                                              outprefix = outprefix,
                                                                              frequency.threshold = 25) else NULL)
    log_message("Done.")
  }
  
  # Creates S3 object containing the updated signature, tables, plots, and survival information
  # Listed parameters are only those directly affecting the scientific results/model
  result <- new_reboot_survival(
    signature = raw_sig_df,
    score_cutoff = score_cutoff,
    expression_data = data,
    clinical_data = clin,
    survival = list(
      univariate = res_logrank,
      multivariate = if (!is.null(multi_cox)) multi_cox else NULL
    ),
    roc_result = if (roc) roc_cutoff else NULL,
    ph_result = if (!is.null(multi_cox)) ph_test else res_logrank$ph,
    plots = plot_list,
    params = list(
      roc = roc,
      multivariate = multivariate,
      variancefilter = variancefilter,
      followup = followup,
      p_cutoff = p.cutoff,
      bootstrap = bootstrap,
      force = force,
      seed = seed
    ),
    metadata = list(
      timestamp = Sys.time(),
      R_version = R.version.string,
      package_version = utils::packageVersion("Reboot")
    ),
    call = match.call()
  )
  
  # Writes execution metadata
  if (saveJSON)
  {
    log_message("Writing JSON metadata...")
    write_reboot_metadata(
      module = "rebootSurvival",
      outprefix = outprefix,
      parameters = list(
        filein = if (is.character(filein)) {filein} else {"data.frame"},
        clinin = if (is.null(clinin)) {NULL} else if (is.character(clinin)) {clinin} else {"data.frame"},
        signature = if (is.character(signature)) {signature} else {class(signature)[1]},
        outprefix = outprefix,
        multivariate = multivariate,
        roc = roc,
        variancefilter = variancefilter,
        followup = followup,
        p.cutoff = p.cutoff,
        bootstrap = bootstrap,
        force = force,
        plots = plots,
        table = table,
        saveJSON = saveJSON,
        saveRDS = saveRDS,
        report = report,
        log = log,
        seed = seed
      ),
      command = "none"
    )
    log_message("Done.")
  }
  
  # Writes RDS object
  if (saveRDS)
  {
    log_message("Writing RDS object...")
    write_reboot_rds(object = result, file = paste0(outprefix, "_survival.rds"))
    log_message("Done.")
  }
  
  # Generates report if requested
  if (!isFALSE(report))
  {
    log_message("Writing REBOOT report...")
    write_reboot_report(object = result, format = report, file = paste0(outprefix, "_report_survival.pdf"))
    log_message("Done.")
  }
  
  # End time counter
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Prints elapsed time
  if (elapsed_time >= 3600) {
    log_message(sprintf("Time to run [module II]: %.2f hours", elapsed_time / 3600))
  } else if (elapsed_time >= 60) {
    log_message(sprintf("Time to run [module II]: %.2f minutes", elapsed_time / 60))
  } else {
    log_message(sprintf("Time to run [module II]: %.2f seconds", elapsed_time))
  }
  
  # Returns the OBJ with the signature and survival results
  if (!log) {log_message("REBOOT survival finished")}
  return(result)
}

#' Run REBOOT complete pipeline
#'
#' @description Runs the complete REBOOT workflow by sequentially executing:
#' \enumerate{
#'   \item \code{rebootRegression()} (module 1)
#'   \item \code{rebootSurvival()} (module 2)
#' }
#'
#' @param filein Character or data.frame. A data.frame or a path to a tab-separated expression/survival file
#' @param outprefix Character. Output prefix used for generated files (Default: \code{"reboot"}).
#' @param bootstrap Integer. Number of bootstraps iterations used in both \code{rebootRegression()} and
#' \code{rebootSurvival()} modules (Default: \code{1}).
#' @param groupsize Integer. Number of features per iteration in \code{rebootRegression()} module (Default: \code{10}).
#' @param percentagefilter Numeric. Correlation filter threshold used in \code{rebootRegression()} (Default: \code{0.3}).
#' @param variancefilter Numeric. Variance filter threshold used in both \code{rebootRegression()}
#' and \code{rebootSurvival()} modules (Default: \code{0.01}).
#' @param followup Numeric. Maximum followup time.
#' @param type Character. Feature type: either \code{"gene"} (default) or \code{"transcript"} used in \code{rebootRegression()}.
#' @param multivariate Logical. If \code{TRUE}, performs multivariate survival analysis
#' in \code{rebootSurvival()} (Default: \code{FALSE}).
#' @param clinin Character or data.frame. A data.frame or a path to a tab-separated clinical file used for multivariate analysis.
#' @param roc Logical. If \code{TRUE}, uses ROC-based cutoff instead of median for scores
#' in \code{rebootSurvival()} (Default: \code{FALSE}).
#' @param p.cutoff Numeric. Univariate covariate filtering cutoff used for multivariate analysis
#' in \code{rebootSurvival()} (Default: \code{0.2}).
#' @param force Logical. If \code{TRUE}, overrides filtering steps related to survival assumptions (Default: \code{FALSE}).
#' @param plots Logical. If \code{TRUE}, generates plots from both modules (Default: \code{FALSE}).
#' @param table Logical. If \code{TRUE}, writes output tables from both modules (Default: \code{FALSE}).
#' @param saveJSON Logical. If \code{TRUE}, saves to disk the metadata in a JSON file (Default: \code{FALSE}).
#' @param saveRDS Logical. If \code{TRUE}, saves to disk the R object in a RDS file (Default: \code{FALSE}).
#' @param report Logical or character. Report generation mode. Use \code{FALSE} (default) to disable report generation,
#' \code{TRUE} or \code{"PDF"} to generate a PDF report, or \code{"HTML"} to generate an HTML report.
#' @param log Logical. If \code{TRUE}, writes log files for both modules (Default: \code{FALSE}).
#' @param seed Integer or \code{NULL}. Random seed controlling every stochastic step in both
#'   \code{rebootRegression()} and \code{rebootSurvival()}, so that results are reproducible across
#'   runs for the same input and parameters (Default: \code{17}). Passed through unchanged to both
#'   modules. Set to \code{NULL} to disable seeding and restore the previous (non-reproducible) behavior.
#'
#' @details
#' The molecular signature generated during the regression step is automatically passed to the survival analysis step.
#'
#' @return An S3 object of class \code{"reboot_complete"} containing:
#' \itemize{
#'   \item \code{regression}: result from \code{rebootRegression()}
#'   \item \code{survival}: result from \code{rebootSurvival()}
#'   \item \code{params}: list of parameters used in the analysis
#'   \item \code{metadata}: metadata information
#'   \item \code{call}: original function call
#' }
#'
#' @seealso
#' Project Github repository: \url{https://github.com/galantelab/reboot} \cr
#'
#' Documentation and tutorials: \url{https://galantelab.github.io/reboot/} \cr
#'
#' Main functions: \code{\link{rebootRegression}} (module I) and \code{\link{rebootSurvival}} (module II)
#'
#' @examples
#' \dontrun{
#' # Example usage
#' result <- rebootComplete(
#'   filein = "expression_data.tsv",
#'   outprefix = "my_analysis",
#'   bootstrap = 100,
#'   groupsize = 10,
#'   percentagefilter = 0.3,
#'   variancefilter = 0.01,
#'   followup = NULL,
#'   type = "gene",
#'   multivariate = TRUE,
#'   clinin = "clinical_data.tsv",
#'   roc = TRUE,
#'   p.cutoff = 0.2,
#'   force = TRUE,
#'   plots = TRUE,
#'   table = TRUE,
#'   saveJSON = TRUE,
#'   saveRDS = TRUE,
#'   report = TRUE,
#'   log = TRUE
#' )
#'
#' # Explore complete analysis results
#' print(result)
#' summary(result)
#' coef(result)
#' plot(result)
#'
#' # Access module-specific results
#' result$regression
#' result$survival
#'
#' # Use module-specific S3 methods
#' summary(result$regression)
#' summary(result$survival)
#' }
#'
#' @export
rebootComplete <- function(filein,
                           outprefix = "reboot",
                           bootstrap = 1,
                           groupsize = 10,
                           percentagefilter = 0.3,
                           variancefilter = 0.01,
                           followup = NULL,
                           type = "gene",
                           multivariate = FALSE,
                           clinin = NULL,
                           roc = FALSE,
                           p.cutoff = 0.2,
                           force = FALSE,
                           plots = FALSE,
                           table = FALSE,
                           saveJSON = FALSE,
                           saveRDS = FALSE,
                           report = FALSE,
                           log = FALSE,
                           seed = 17)
{
  # Start time counter
  start_time <- Sys.time()

  # Checks SEED argument (module I and II each set/restore it locally around their own work)
  .reboot_check_seed(seed)
  if (!log) {log_message("REBOOT complete log started")}
  
  # Runs module 1 (regression)
  log_message("=============== REBOOT: module I - regression ===============")
  regression_result <- rebootRegression(
    filein = filein,
    outprefix = outprefix,
    bootstrap = bootstrap,
    groupsize = groupsize,
    percentagefilter = percentagefilter,
    variancefilter = variancefilter,
    followup = followup,
    type = type,
    force = force,
    plots = plots,
    table = table,
    saveJSON = FALSE,
    saveRDS = FALSE,
    report = FALSE,
    log = log,
    seed = seed
  )
  log_message("============================================================")
  
  # Runs module 2 (survival)
  log_message("=============== REBOOT: module II - survival ===============")
  survival_result <- rebootSurvival(
    filein = filein,
    signature = regression_result,
    outprefix = outprefix,
    multivariate = multivariate,
    clinin = clinin,
    roc = roc,
    variancefilter = variancefilter,
    followup = followup,
    p.cutoff = p.cutoff,
    bootstrap = bootstrap,
    force = force,
    plots = plots,
    table = table,
    saveJSON = FALSE,
    saveRDS = FALSE,
    report = FALSE,
    log = log,
    seed = seed
  )
  log_message("============================================================")
  
  # Creates final S3 object
  # Listed parameters are only those directly affecting the scientific results/model
  result <- new_reboot_complete(
    regression = regression_result,
    survival = survival_result,
    params = list(
      bootstrap = bootstrap,
      groupsize = groupsize,
      percentagefilter = percentagefilter,
      variancefilter = variancefilter,
      followup = followup,
      type = type,
      multivariate = multivariate,
      roc = roc,
      p_cutoff = p.cutoff,
      force = force,
      seed = seed
    ),
    metadata = list(
      timestamp = Sys.time(),
      R_version = R.version.string,
      package_version = utils::packageVersion("Reboot")
    ),
    call = match.call()
  )
  
  # Writes execution metadata
  if (saveJSON)
  {
    log_message("Writing JSON metadata...")
    write_reboot_metadata(
      module = "rebootComplete",
      outprefix = outprefix,
      parameters = list(
        filein = if (is.character(filein)) {filein} else {"data.frame"},
        clinin = if (is.null(clinin)) {NULL} else if (is.character(clinin)) {clinin} else {"data.frame"},
        signature = "data.frame",
        outprefix = outprefix,
        multivariate = multivariate,
        roc = roc,
        groupsize = groupsize,
        percentagefilter = percentagefilter,
        variancefilter = variancefilter,
        followup = followup,
        p.cutoff = p.cutoff,
        bootstrap = bootstrap,
        type = type,
        force = force,
        plots = plots,
        table = table,
        saveJSON = saveJSON,
        saveRDS = saveRDS,
        report = report,
        log = log,
        seed = seed
      ),
      command = "none"
    )
    log_message("Done.")
  }
  
  # Writes RDS object
  if (saveRDS)
  {
    log_message("Writing RDS object...")
    write_reboot_rds(object = result, file = paste0(outprefix, "_complete.rds"))
    log_message("Done.")
  }
  
  # Generates report if requested
  if (!isFALSE(report))
  {
    log_message("Writing REBOOT report...")
    write_reboot_report(object = result, format = report, file = paste0(outprefix, "_report_complete.pdf"))
    log_message("Done.")
  }
  
  # End time counter
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Prints elapsed time
  if (elapsed_time >= 3600) {
    log_message(sprintf("Time to run [complete workflow]: %.2f hours", elapsed_time / 3600))
  } else if (elapsed_time >= 60) {
    log_message(sprintf("Time to run [complete workflow]: %.2f minutes", elapsed_time / 60))
  } else {
    log_message(sprintf("Time to run [complete workflow]: %.2f seconds", elapsed_time))
  }
  
  # Returns object
  return(result)
}
