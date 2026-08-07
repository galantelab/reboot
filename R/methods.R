###### REGRESSION (MODULE I) ######

#' Print of REBOOT signature
#'
#' @description Displays a 'print' of the `reboot_signature` object, including the number of
#' selected features and a preview of the top entries.
#'
#' @param x An object of class \code{reboot_signature}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#'
#' @method print reboot_signature
#' @export
print.reboot_signature <- function(x, ...) {
  
  # Checks if input object is of correct type
  if (!inherits(x, "reboot_signature")) {stop("Object must be of class 'reboot_signature'")}
  
  cat("REBOOT signature object\n")
  cat("Number of features:", nrow(x$signature), "\n\n")
  
  # Checks if signature is not empty
  if (nrow(x$signature) == 0)
  {
    cat("No features selected.\n")
    return(invisible(x))
  }
  
  # Prints top 10 features
  utils::head(x$signature, 10)
  if (nrow(x$signature) > 10) {cat("\n... (truncated)\n")}
  invisible(x)
}

#' Summary of REBOOT signature
#'
#' @description Computes summary statistics for a `reboot_signature` object, including the
#' number of features and descriptive statistics of coefficients.
#'
#' @param object An object of class \code{reboot_signature}.
#' @param ... Additional arguments (not used).
#'
#' @return A list of class `summary.reboot_signature` containing:
#' \describe{
#'   \item{n_features}{Number of features in the signature}
#'   \item{coef_min}{Minimum coefficient value}
#'   \item{coef_max}{Maximum coefficient value}
#'   \item{coef_median}{Median coefficient value}
#'   \item{coef_mean}{Mean coefficient value}
#' }
#'
#'
#' @method summary reboot_signature
#' @export
summary.reboot_signature <- function(object, ...) {
  
  # Checks if input object is of correct type
  if (!inherits(object, "reboot_signature")) {stop("Object must be of class 'reboot_signature'")}
  
  # Extracts only the coefficients
  coefs <- object$signature$coefficient
  
  # Checks if signature is not empty and prints summary statistics
  if (length(coefs) == 0) {
    out <- list(
      n_features = 0,
      coef_min = NA_real_,
      coef_max = NA_real_,
      coef_median = NA_real_,
      coef_mean = NA_real_
    )
  } else {
    out <- list(
      n_features = length(coefs),
      coef_min = base::min(coefs, na.rm = TRUE),
      coef_max = base::max(coefs, na.rm = TRUE),
      coef_median = stats::median(coefs, na.rm = TRUE),
      coef_mean = base::mean(coefs, na.rm = TRUE)
    )
  }
  
  class(out) <- "summary.reboot_signature"
  return(out)
}

#' Print summary of REBOOT signature
#'
#' @description Custom print method for objects of class `summary.reboot_signature`.
#' Displays key statistics of the signature in a user-friendly format.
#'
#' @param x An object of class \code{summary.reboot_signature}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.summary.reboot_signature <- function(x, ...) {
  
  # Checks if input object is of correct type
  if (!inherits(x, "summary.reboot_signature")) {stop("Object must be of class 'summary.reboot_signature'")}
  
  cat("REBOOT signature summary\n")
  cat("Number of features:", x$n_features, "\n")
  
  # Checks if signature summary statistics are not NA and prints them
  if (is.na(x$coef_min) || is.na(x$coef_max)) {
    cat("Coefficient range: NA\n")
  } else {
    cat("Coefficient range:", x$coef_min, "to", x$coef_max, "\n")
  }
  cat("Median coefficient:", x$coef_median, "\n")
  cat("Mean coefficient:", x$coef_mean, "\n")
  invisible(x)
}

#' Convert REBOOT signature to data.frame
#'
#' @description Coerces a `reboot_signature` object into a standard `data.frame`,
#' returning the underlying signature table.
#'
#' @param x An object of class \code{reboot_signature}.
#' @param ... Additional arguments (not used).
#'
#' @return A `data.frame` containing the signature with 'features' and 'coefficients'.
#'
#' @method as.data.frame reboot_signature
#' @export
as.data.frame.reboot_signature <- function(x, ...) {
  
  # Checks if input object is of correct type
  if (!inherits(x, c("reboot_signature", "reboot_survival")))
  {
    stop("Object must be of class 'reboot_signature' or 'reboot_survival'")
  }
  
  # Checks if 'reboot_signature' is valid
  if (is.null(x$signature) || !is.data.frame(x$signature)) {
    stop("Invalid 'reboot' class object: missing signature data")
  }
  
  # Converts signature to 'data.frame'
  df <- x$signature
  rownames(df) <- NULL
  return(df)
}

#' Plot REBOOT signature
#'
#' @description Generates visualization(s) for a `reboot_signature` object.
#'
#' @param x An object of class \code{reboot_signature}.
#' @param type Character. Type of plot: \code{"lolli"}, \code{"hist"}, or \code{"all"} (default).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A \code{ggplot} object or a named list of plots.
#'
#' @method plot reboot_signature
#' @export
plot.reboot_signature <- function(x, type = "all", ...)
{
  # Checks if input object is of correct type
  if (!inherits(x, "reboot_signature")) {stop("Object must be of class 'reboot_signature'")}
  
  # Validates type argument
  if (!is.character(type) || length(type) != 1) {stop("'type' must be a single character string")}
  
  # Extracts signature
  sig <- x$signature
  
  # Lollipop plot
  if (type == "lolli") {return(reboot_lollipop(signature = sig, ...))}
  
  # Histogram
  if (type == "hist") {return(reboot_histogram(signature = sig, ...))}
  
  # All plots
  if (type == "all")
  {
    plot_list <- list()
    
    plot_list$lollipop <- reboot_lollipop(signature = sig, ...)
    print(plot_list$lollipop)
    
    plot_list$histogram <- reboot_histogram(signature = sig, ...)
    print(plot_list$histogram)
    
    return(invisible(plot_list))
  }
  
  stop("Invalid 'type'. Must be 'lolli', 'hist', or 'all' (default)")
}

#' Coefficients from REBOOT signature
#'
#' @description Returns the coefficients from a `reboot_signature` object.
#'
#' @param object An object of class \code{reboot_signature}.
#' @param ... Additional arguments (not used).
#'
#' @return A named numeric vector of coefficients.
#'
#' @method coef reboot_signature
#' @export
coef.reboot_signature <- function(object, ...) {
  
  # Checks if input object is of correct type
  if (!inherits(object, "reboot_signature")) {stop("Object must be of class 'reboot_signature'")}
  
  # Gets the signature
  sig <- object$signature
  
  # Checks if 'reboot_signature' is valid
  if (is.null(sig) || !all(c("feature", "coefficient") %in% colnames(sig))) {
    stop("Invalid 'reboot_signature' object: missing signature data")
  }
  
  # Creates named vector of coefficients similar to what's seem in coef(glm) or coef(lm)
  valid <- !is.na(sig$coefficient)
  coefs <- sig$coefficient[valid]
  names(coefs) <- sig$feature[valid]
  
  return(coefs)
}

###### SURVIVAL (MODULE II) ######

#' Convert REBOOT survival object to data.frame
#'
#' @description Coerces a `reboot_survival` object into a standard `data.frame`, returning the underlying signature table.
#'
#' @param x An object of class \code{reboot_survival}.
#' @param ... Additional arguments (not used).
#'
#' @return A `data.frame` containing the signature with 'features' and 'coefficients'.
#'
#' @method as.data.frame reboot_survival
#' @export
as.data.frame.reboot_survival <- as.data.frame.reboot_signature

#' Print REBOOT survival object
#'
#' @description Displays a concise summary of a `reboot_survival` object, including signature, score, and survival analyses.
#'
#' @param x An object of class \code{reboot_survival}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#'
#' @method print reboot_survival
#' @export
print.reboot_survival <- function(x, ...) {
  
  # Checks if input object is of correct type
  if (!inherits(x, "reboot_survival")) {stop("Object must be of class 'reboot_survival'")}
  
  cat("REBOOT survival object\n")
  cat("-----------------------\n")
  
  # Signature
  cat("Signature features:", nrow(x$signature), "\n")
  
  # Samples
  if (!is.null(x$expression_data)) {cat("Samples:", nrow(x$expression_data), "\n")}
  
  # Cutoff
  if (!is.null(x$score_cutoff)) {cat("Score cutoff:", signif(x$score_cutoff, 4), "\n")}
  
  # ROC
  cat("ROC analysis:", ifelse(is.null(x$roc_result), "No", "Yes"), "\n")
  
  # Multivariate
  has_multi <- !is.null(x$survival$multivariate)
  cat("Multivariate analysis:", ifelse(has_multi, "Yes", "No"), "\n")
  
  # Clinical data
  cat("Clinical data:", ifelse(is.null(x$clinical_data), "No", "Yes"), "\n")
  
  invisible(x)
}

#' Summary of REBOOT survival object
#'
#' @description Computes summary statistics and analytical results from a `reboot_survival` object.
#'
#' @param object An object of class \code{reboot_survival}.
#' @param ... Additional arguments (not used).
#'
#' @return A list of class \code{summary.reboot_survival}.
#'
#' @method summary reboot_survival
#' @export
summary.reboot_survival <- function(object, ...) {
  
  # Checks if input object is of correct type
  if (!inherits(object, "reboot_survival")) {stop("Object must be of class 'reboot_survival'")}
  
  # Basic information
  out <- list(
    n_features = nrow(object$signature),
    n_samples = if (!is.null(object$expression_data)) {nrow(object$expression_data)} else {NA_integer_},
    cutoff = object$score_cutoff,
    score_summary = if (!is.null(object$expression_data$score)) {summary(object$expression_data$score)} else {NULL},
    score_groups = if (!is.null(object$clinical_data$score_group))
    {
      table(object$clinical_data$score_group)
    } else {NULL},
    roc_auc = if (!is.null(object$roc_result$auc)) {object$roc_result$auc} else {NA_real_},
    univariate_pvalue = if (!is.null(object$survival$univariate$table))
    {
      object$survival$univariate$table$log.rank.pvalue
    } else {NA_real_},
    multivariate = !is.null(object$survival$multivariate),
    n_clinical_variables = if (!is.null(object$clinical_data)) {ncol(object$clinical_data) - 4} else {0}
  )
  
  class(out) <- "summary.reboot_survival"
  
  return(out)
}

#' Print summary of REBOOT survival object
#'
#' @description Custom print method for objects of class `summary.reboot_survival`.
#'
#' @param x An object of class \code{summary.reboot_survival}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.summary.reboot_survival <- function(x, ...) {
  
  # Checks if input object is of correct type
  if (!inherits(x, "summary.reboot_survival")) {stop("Object must be of class 'summary.reboot_survival'")}
  
  cat("REBOOT survival summary\n")
  cat("------------------------\n")
  
  cat("Features:", x$n_features, "\n")
  cat("Samples:", x$n_samples, "\n")
  cat("Score cutoff:", signif(x$cutoff, 4), "\n\n")
  
  cat("Score distribution:\n")
  print(x$score_summary)
  
  cat("\nScore groups:\n")
  print(x$score_groups)
  
  cat("\nUnivariate log-rank p-value:", signif(x$univariate_pvalue, 4), "\n")
  
  if (!is.na(x$roc_auc)) {cat("ROC AUC:", signif(x$roc_auc, 4), "\n")}
  
  cat("Multivariate analysis:", ifelse(x$multivariate, "Yes", "No"), "\n")
  
  cat("Clinical covariates:", x$n_clinical_variables, "\n")
  
  invisible(x)
}

#' Coefficients from REBOOT survival object
#'
#' @description Returns the coefficients from the molecular signature stored inside a `reboot_survival` object.
#'
#' @param object An object of class \code{reboot_survival}.
#' @param ... Additional arguments (not used).
#'
#' @return A named numeric vector of coefficients.
#'
#' @method coef reboot_survival
#' @export
coef.reboot_survival <- function(object, ...) {
  
  # Checks if input object is of correct type
  if (!inherits(object, "reboot_survival")) {stop("Object must be of class 'reboot_survival'")}
  
  sig <- object$signature
  
  # Checks signature structure
  if (is.null(sig) || !all(c("feature", "coefficient") %in% colnames(sig)))
  {
    stop("Invalid 'reboot_survival' object: missing signature data")
  }
  
  # Creates named coefficient vector
  valid <- !is.na(sig$coefficient)
  coefs <- sig$coefficient[valid]
  names(coefs) <- sig$feature[valid]
  
  return(coefs)
}

#' Plot REBOOT survival object
#'
#' @description Generates visualization(s) from a `reboot_survival` object.
#'
#' @param x An object of class \code{reboot_survival}.
#' @param type Character. Plot type:
#' \code{"km"}, \code{"roc"}, \code{"forest"},
#' \code{"ph"}, \code{"bar"}, or \code{"all"} (default).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A \code{ggplot} object or a named list of plots.
#'
#' @method plot reboot_survival
#' @export
plot.reboot_survival <- function(x, type = "all", ...)
{
  # Checks object type
  if (!inherits(x, "reboot_survival")) {stop("Object must be of class 'reboot_survival'")}
  
  # Validates type argument
  if (!is.character(type) || length(type) != 1) {stop("'type' must be a single character string")}
  
  # Defines initial available plot types
  valid_types <- c("all", "km", "ph")
  
  # Optional ROC
  if (!is.null(x$roc_result)) {valid_types <- c(valid_types, "roc")}
  
  # Optional multivariate plots
  if (!is.null(x$survival$multivariate))
  {
    valid_types <- c(valid_types, "forest")
    if (!is.null(x$survival$multivariate$bootstrap.table)) {valid_types <- c(valid_types, "bar")}
  }
  
  # Validates requested plot type
  if (!type %in% valid_types) {stop("Invalid 'type'. Must be one of: ", paste(valid_types, collapse = ", "))}
  
  # KM plot
  if (type == "km")
  {
    return(reboot_km(data = x$expression_data, cutoff = x$score_cutoff,
                     pval = x$survival$univariate$table$log.rank.pvalue, ...))
  }
  
  # PH plot
  if (type == "ph") {return(reboot_phplot(x = x$ph_result, is_multi = !is.null(x$survival$multivariate), ...))}
  
  # ROC plot
  if (type == "roc") {return(reboot_roc(x = x$roc_result, ...))}
  
  # Forest plot
  if (type == "forest") {return(reboot_forest(model_object = x$survival$multivariate$model, ...))}
  
  # Bootstrap barplot
  if (type == "bar") {return(reboot_barplot(plot_df = x$survival$multivariate$bootstrap.table, ...))}
  
  # All plots using controlled recursion
  if (type == "all")
  {
    plot_s3_list <- list()
    plot_s3_list$km <- plot(x, type = "km", ...)
    print(plot_s3_list$km)
    plot_s3_list$ph <- plot(x, type = "ph", ...)
    print(plot_s3_list$ph)
    if ("roc" %in% valid_types)
    {
      plot_s3_list$roc <- plot(x, type = "roc", ...)
      grDevices::replayPlot(plot_s3_list$roc)
    }
    if ("forest" %in% valid_types)
    {
      plot_s3_list$forest <- plot(x, type = "forest", ...)
      print(plot_s3_list$forest)
    }
    if ("bar" %in% valid_types)
    {
      plot_s3_list$bar <- plot(x, type = "bar", ...)
      print(plot_s3_list$bar)
    }
    return(invisible(plot_s3_list))
  }
}

###### COMPLETE (MODULES I + II) ######

#' Print a REBOOT complete analysis object
#'
#' @description Prints a concise summary of a `reboot_complete` object,
#' including information from both regression and survival modules.
#'
#' @param x Object of class \code{"reboot_complete"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @method print reboot_complete
#' @export
print.reboot_complete <- function(x, ...)
{
  cat("\n")
  cat("===== REBOOT Complete Analysis =====\n\n")
  
  # Regression module
  cat("[Module 1 - Regression]\n")
  cat("Selected features :", nrow(x$regression$signature), "\n")
  cat("Feature type :", x$regression$params$type, "\n")
  cat("Bootstraps :", x$regression$params$bootstrap, "\n\n")
  
  # Survival module
  cat("[Module 2 - Survival]\n")
  cat("Samples :", nrow(x$survival$expression_data), "\n")
  cat("ROC :", x$survival$params$roc, "\n")
  cat("Multivariate :", x$survival$params$multivariate, "\n")
  
  # Clinical covariates - removing survival/score columns
  if (!is.null(x$survival$clinical_data))
  {
    n_covars <- ncol(x$survival$clinical_data) - 4
    cat("Clinical variables:", n_covars, "\n")
  }
  
  # Multivariate status
  multi_status <- !is.null(x$survival$survival$multivariate)
  cat("Multi-Cox model :", ifelse(multi_status, "available", "not available"), "\n")
  
  # Metadata
  cat("\n")
  cat("Generated :", as.character(x$metadata$timestamp), "\n")
  cat("R version :", x$metadata$R_version, "\n")
  cat("Package :", paste0("Reboot v", x$metadata$package_version), "\n")
  
  invisible(x)
}

#' Summarize a REBOOT complete analysis object
#'
#' @description Generates a detailed summary of a `reboot_complete` object,
#' integrating outputs from both regression and survival modules.
#'
#' @param object Object of class \code{"reboot_complete"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class \code{"summary.reboot_complete"}.
#'
#' @method summary reboot_complete
#' @export
summary.reboot_complete <- function(object, ...)
{
  # Regression summary
  regression_summary <- summary(object$regression)
  
  # Survival summary
  survival_summary <- summary(object$survival)
  
  # Creates summary object
  result <- list(
    regression = regression_summary,
    survival = survival_summary,
    params = object$params,
    metadata = object$metadata,
    call = object$call
  )
  
  class(result) <- "summary.reboot_complete"
  
  return(result)
}

#' Print summary of a REBOOT complete analysis object
#'
#' @description Prints the summary generated by \code{summary.reboot_complete()}.
#'
#' @param x Object of class \code{"summary.reboot_complete"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.summary.reboot_complete <- function(x, ...)
{
  cat("\n")
  cat("===== REBOOT Complete Summary =====\n")
  
  cat("\n--- Regression module ---\n")
  print(x$regression)
  
  cat("\n--- Survival module ---\n")
  print(x$survival)
  
  cat("\n--- General parameters ---\n")
  print(x$params)
  
  cat("\n--- Metadata ---\n")
  print(x$metadata)
  
  invisible(x)
}

#' Extract coefficients from a REBOOT complete analysis object
#'
#' @description Extracts feature coefficients from the survival module of a `reboot_complete` object.
#'
#' @param object Object of class \code{"reboot_complete"}.
#' @param ... Additional arguments passed to downstream methods.
#'
#' @return A named numeric vector of coefficients.
#'
#' @method coef reboot_complete
#' @export
coef.reboot_complete <- function(object, ...)
{
  stats::coef(object$survival, ...)
}

#' Plot a REBOOT complete analysis object
#'
#' @description Displays plots generated during the complete REBOOT workflow.
#'
#' @param x Object of class \code{"reboot_complete"}.
#' @param type Character. Vector indicating which module plots to display. Possible values are:
#' \code{"regression"}, \code{"survival"}, or \code{"all"} (default).
#' @param ... Additional arguments passed to module-specific plot methods.
#'
#' @return Invisibly returns a list containing the generated plots
#'
#' @method plot reboot_complete
#' @export
plot.reboot_complete <- function(x, type = "all", ...)
{
  valid <- c("all", "regression", "survival")
  if (!all(type %in% valid)) {stop("Argument 'type' must be: 'all', 'regression', or 'survival'")}
  
  plot_list <- list()
  
  # Regression module
  if ("all" %in% type || "regression" %in% type) {plot_list$regression <- plot(x$regression, ...)}
  
  # Survival module
  if ("all" %in% type || "survival" %in% type) {plot_list$survival <- plot(x$survival, ...)}
  
  invisible(plot_list)
}
