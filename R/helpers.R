utils::globalVariables(c("variable"))

###### MAIN (REGRESSION - MODULE I) ######

#' Applies maximum follow-up censoring to survival data
#'
#' @description Applies administrative censoring to survival data by limiting follow-up time to a user-defined
#' maximum value. Samples with survival times greater than the specified follow-up are censored.
#'
#' @param data Data.frame containing survival information.
#' @param followup Numeric. Maximum follow-up time. If \code{NULL} (default), no filtering is applied.
#'
#' @return A filtered data.frame with adjusted \code{OS} and \code{OS.time}.
#' @keywords internal
reboot_followup <- function(data, followup = NULL) {
  # Validates input data
  if (!is.data.frame(data)) {stop("'data' must be a data.frame")}
  required_cols <- c("OS", "OS.time")
  if (!all(required_cols %in% colnames(data))) {stop("Input data must contain columns: 'OS' and 'OS.time'")}
  
  # Returns original data if followup is NULL
  if (is.null(followup))
  {
    message("Returning original data.frame (no followup max value passed)...")
    return(data)
  }
  
  # Validates followup
  if (!is.numeric(followup) || length(followup) != 1 || followup <= 0) {stop("'followup' must be a positive numeric value")}
  
  # Skips filtering if followup exceeds maximum observed time
  max_followup <- max(data$OS.time, na.rm = TRUE)
  if (followup >= max_followup)
  {
    message("Returning original data.frame (followup is greater than or equal to the max 'OS.time')...")
    return(data)
  }
  
  # Identifies samples exceeding followup threshold
  idx <- which(data$OS.time > followup)
  
  # Applies administrative censoring
  data$OS[idx] <- 0
  data$OS.time[idx] <- followup
  
  return(data)
}

#' Variance filtering for REBOOT pipeline
#'
#' @description Internal helper function used by \code{rebootRegression()}.
#' Applies variance-based filtering to remove low-variance features after normalization.
#' Also performs basic checks on survival data.
#'
#' @param cmatrix A data.frame with survival data (first two columns) and expression values.
#' @param var_threshold Numeric. Variance threshold for filtering.
#' @param force Logical. Whether to bypass survival-related checks.
#'
#' @return A filtered data.frame.
#' @keywords internal
reboot_varfun <- function(cmatrix, var_threshold, force) {
  
  # Handle missing values
  if (anyNA(cmatrix)) {
    message("Checking NAs and imputing with mice...")
    impu <- mice::mice(cmatrix, print = FALSE)
    cmatrix <- mice::complete(impu)
  }
  
  # Forces the first 2 columns to be named OS and OS.time
  colnames(cmatrix)[1:2] <- c("OS", "OS.time")
  
  # Handles cases where the user has only 1 feature
  if (ncol(cmatrix) == 3) {
    if (any(cmatrix[[3]] == 0)) {stop("Feature column contains zero values, cannot normalize")}
    divisor <- max(cmatrix[[3]])
    dividendo <- cmatrix[[3]]
    normalized <- divisor / dividendo
    variance <- stats::var(normalized)
    
    if (variance < var_threshold) {stop("All columns rejected by variance filter")}
    message("No columns rejected by variance filter")
  } else {
    maxes <- apply(cmatrix[, 3:ncol(cmatrix), drop = FALSE], 2, max)
    
    if (0 %in% maxes) {stop("Columns containing only zeros found in input data. Remove such columns and try again.")}
    message("Calculating normalized variances")
    
    divisor <- cmatrix[, 3:ncol(cmatrix), drop = FALSE]
    dividendo <- matrix(maxes, nrow = nrow(cmatrix), ncol = length(maxes), byrow = TRUE)
    normalized <- divisor / dividendo
    variances <- apply(normalized, 2, stats::var)
    filtered <- which(variances > var_threshold)
    
    if (length(filtered) == 0) {
      stop("All columns rejected by variance filter")
    } else if (length(filtered) == length(variances)) {
      message("No columns rejected by variance filter")
    } else {
      losers <- names(variances)[-filtered]
      cmatrix <- cmatrix[, c(1, 2, filtered + 2), drop = FALSE]
      message(length(losers), " columns removed by variance filter: ", paste(losers, collapse = ", "))
    }
    
    # Survival checks
    if (!force) {
      OSstatus <- cmatrix[, 1]
      percentage <- sum(OSstatus) / length(OSstatus)
      
      if (percentage < 0.2 || percentage > 0.8) {stop("Survival status proportion is not adequate for analysis")}
      
      followup <- cmatrix[, 2]
      uplimit <- max(followup)
      normalized <- followup / uplimit
      fvar <- stats::var(normalized)
      
      if (fvar < var_threshold) {stop("Follow-up variance did not pass the threshold")}
    }
  }
  
  return(cmatrix)
}

#' Schoenfeld test filtering for REBOOT pipeline
#'
#' @description Internal helper function used by \code{rebootRegression()}.
#' Filters features that violate proportional hazards assumption using Schoenfeld residuals.
#'
#' @param full_data A data.frame with survival data and features.
#'
#' @return A filtered data.frame.
#' @keywords internal
reboot_ph_assumptions <- function(full_data) {
  
  filt <- character()
  attributes <- colnames(full_data)[-c(1, 2)]
  
  for (i in attributes) {
    formula_str <- paste("survival::Surv(OS.time, OS) ~", i)
    
    # Modification to run Cox model faster (avoid redundancy)
    tryCatch({
      phmodel <- survival::coxph(stats::as.formula(formula_str), data = full_data)
      
      if (!is.na(phmodel$coef)) {
        schoen <- survival::cox.zph(phmodel)
        pval <- schoen$table[1, 3]
        if (!is.na(pval) && pval > 0.05) {filt <- c(filt, i)}
      }
    },
    warning = function(w) {message("Warning in feature ", i, ": ", conditionMessage(w))},
    error = function(e) {message("Skipping feature ", i, " due to error: ", conditionMessage(e))})
  }
  
  losers <- setdiff(attributes, filt)
  
  # Controls user output to avoid pollution
  if (length(losers) > 50) {
    message(length(losers), " columns removed by Schoenfeld test (showing first 50): ",
            paste(utils::head(losers, 50), collapse = ", "))
  } else {
    message(length(losers), " columns removed by Schoenfeld test: ", paste(losers, collapse = ", "))
  }
  
  if (length(filt) == 0) {stop("No features passed the Schoenfeld test")}
  
  cols <- match(filt, colnames(full_data))
  cols <- cols[!is.na(cols)]
  return(full_data[, c(1, 2, cols), drop = FALSE])
}

#' Feature count validation for REBOOT pipeline
#' Substitutes the old numberfilter1() and numberfilter2() functions
#'
#' @description Internal helper function used by \code{rebootRegression()}.
#' Checks if sufficient features remain for bootstrap regression.
#' If not, performs a single regression and returns a signature.
#'
#' @param dataf A data.frame with survival data and features.
#' @param g Integer. Group size.
#'
#' @return NULL if pipeline should continue, or a data.frame signature.
#' @keywords internal
reboot_feature_check <- function(dataf, g) {
  
  n_features <- ncol(dataf) - 2
  
  # Extreme case: no features
  if (n_features <= 0) {stop("No features remaining after filtering steps")}
  
  # Low features case: direct regression
  if (n_features < g) {
    message("Number of features is smaller than group size. Performing single multivariate regression...")
    coefs <- reboot_regression(dataf)
    if (is.null(coefs) || length(coefs) == 0 || all(coefs == 0)) {stop("No signature found: all coefficients are zero")}
    signature <- data.frame(feature = names(coefs), coefficient = unname(coefs), stringsAsFactors = FALSE)
    signature$feature <- gsub("__", "-", signature$feature)
    return(signature)
  }
  # Normal expected case: move forward in the pipeline
  return(NULL)
}

#' Bootstrap regression for REBOOT pipeline
#'
#' @description Internal helper function used by \code{rebootRegression()}.
#' Performs bootstrap resampling, applies correlation filtering and regression,
#' and aggregates coefficients across iterations to derive a final signature.
#'
#' @param full_data A data.frame with survival data and features.
#' @param n_boot Integer. Number of bootstrap iterations.
#' @param group_size Integer. Number of features sampled per iteration.
#' @param cor_threshold Numeric. Correlation filter threshold.
#' @param coef_threshold Numeric. Minimum absolute coefficient value to retain features.
#'
#' @return A data.frame with the following signature columns:
#' \describe{
#'   \item{feature}{Feature name (gene/transcript)}
#'   \item{coefficient}{Mean coefficient across bootstrap iterations}
#'   \item{sd}{Standard deviation of coefficients}
#' }
#' @keywords internal
reboot_bootstrapfun <- function(full_data, n_boot, group_size, cor_threshold, coef_threshold) {
  
  message("Starting bootstrap with ", n_boot, " iterations...")
  features <- colnames(full_data)[-c(1, 2)]
  results <- stats::setNames(vector("list", length(features)), features)
  i <- 1
  
  # Avoids infinite loop
  attempts <- 0
  max_attempts <- n_boot * 10
  
  while (i <= n_boot) {
    attempts <- attempts + 1
    if (attempts > max_attempts) {stop("Too many failed bootstrap attempts (correlation filter too strict?)")}
    
    # Avoids spam messages
    if (i %% 10 == 0) {message("Processing iteration: ", i)}
    
    # Subsampling
    if (ncol(full_data) == 3) {
      cmatrix <- full_data
    } else {
      # NOTE: no per-iteration 'seed' is passed here on purpose. Reproducibility across the whole
      # bootstrap loop (and the LASSO cross-validation folds and imputation steps that follow) is
      # controlled once, upstream, by rebootRegression()'s 'seed' argument. Reseeding to a small
      # sequential integer on every iteration would make draws less random and would not cover the
      # other stochastic steps in the pipeline (see rebootRegression() for details).
      cmatrix <- reboot_subsample(full_data, group_size)
      if (group_size > 1 && reboot_corfun(cmatrix, cor_threshold) == 1) {next} # Correlation filter
    }
    
    # Regression
    coefs <- reboot_regcall(cmatrix, group_size, full_data)
    
    if (!is.null(coefs)) {
      if (is.null(names(coefs))) {stop("reboot_regcall() must return a named vector")}
      coef_names <- names(coefs)
      for (j in seq_along(coefs)) {
        fname <- coef_names[j]
        results[[fname]] <- c(results[[fname]], coefs[[j]])
      }
      i <- i + 1
    }
  }
  
  # Aggregates results
  summary_list <- lapply(names(results), function(feature) {
    vals <- results[[feature]]
    if (length(vals) == 0) return(NULL)
    data.frame(feature = feature, coefficient = mean(vals, na.rm = TRUE),
               sd = if (sum(!is.na(vals)) > 1) stats::sd(vals, na.rm = TRUE) else NA_real_,
               stringsAsFactors = FALSE)
  })
  
  summary_list <- Filter(Negate(is.null), summary_list)
  tt <- do.call(rbind, summary_list)
  
  # Checks validity of results
  if (is.null(tt) || nrow(tt) == 0) {stop("No coefficients estimated during bootstrap")}
  if (anyNA(tt$coefficient)) {warning("NA coefficients detected. Consider increasing bootstrap iterations.")}
  
  # Applies coefficient threshold
  tt <- tt[abs(tt$coefficient) >= coef_threshold, , drop = FALSE]
  if (nrow(tt) == 0 || all(tt$coefficient == 0)) {stop("No signature found: coefficients below threshold")}
  tt$feature <- gsub("__", "-", tt$feature)
  message("Bootstrap completed successfully.")

  return(tt)
}

###### SUPPORT (REGRESSION - MODULE I) ######

#' Correlation filtering for REBOOT pipeline
#'
#' @description Internal helper function used by \code{rebootRegression()}.
#' Evaluates pairwise Spearman correlation between features and determines
#' whether the current bootstrap sample should be discarded based on a threshold.
#'
#' @param cmatrix A data.frame with survival data and features.
#' @param cor_threshold Numeric. Maximum allowed proportion of highly correlated feature pairs.
#'
#' @return Integer. Returns 1 if the sample should be discarded, 0 otherwise.
#' @keywords internal
reboot_corfun <- function(cmatrix, cor_threshold) {
  
  n_features <- ncol(cmatrix) - 2
  
  # No correlation possible with < 2 features
  if (n_features < 2) return(0)
  
  # Optimized (and faster) way of creating vectors (pre-defining types and length)
  feature_idx <- 3:ncol(cmatrix)
  n_pairs <- choose(length(feature_idx), 2)
  correlations <- numeric(n_pairs)
  pvals <- numeric(n_pairs)
  pair_names <- character(n_pairs)
  k <- 1
  
  # Performs correlations
  for (i in seq_len(length(feature_idx) - 1)) {
    for (j in (i + 1):length(feature_idx)) {
      t <- feature_idx[i]
      u <- feature_idx[j]
      
      # Brings more security to correlation
      test <- tryCatch(
        stats::cor.test(cmatrix[, t], cmatrix[, u], method = "spearman", exact = FALSE),
        error = function(e) return(NULL)
      )
      if (is.null(test)) next
      
      # No rounding here to make more precise decisions
      correlations[k] <- as.numeric(test$estimate)
      pvals[k] <- test$p.value
      pair_names[k] <- paste(colnames(cmatrix)[t], colnames(cmatrix)[u], sep = "_")
      
      k <- k + 1
    }
  }
  
  # Avoids division by zero
  if (length(pvals) == 0) return(0)
  
  high_corr <- (correlations > 0.80) & (pvals < 0.05)
  proportion <- sum(high_corr) / length(pvals)
  
  if (proportion >= cor_threshold) {
    message("Iteration skipped due to high correlation among features.")
    if (any(high_corr)) {message("High correlation pairs: ", paste(pair_names[high_corr], collapse = ", "))}
    return(1)
  }
  
  return(0)
}

#' Bootstrap subsampling for REBOOT pipeline
#'
#' @description Internal helper function used by \code{rebootRegression()}.
#' Randomly samples a subset of features for each bootstrap iteration,
#' preserving survival columns and enabling stochastic feature selection.
#'
#' @param full_data A data.frame with survival data and features
#' @param group_size Integer. Number of features to sample in the current iteration.
#' @param seed Integer or NULL. Optional random seed for reproducibility. If NULL, no seed is set.
#'
#' @return A data.frame containing the survival data and the sampled features.
#' @keywords internal
reboot_subsample <- function(full_data, group_size, seed = NULL) {
  
  # Checks validity of input
  if (ncol(full_data) <= 2) {stop("No features available for subsampling")}
  if (!is.null(seed)) {set.seed(seed)}
  feature_cols <- colnames(full_data)[3:ncol(full_data)]
  if (group_size > length(feature_cols)) {stop("'group_size' is larger than available features")}
  
  # Makes bootstrap resampling with no replacement
  shuffle <- sample(feature_cols, size = group_size, replace = FALSE)
  message("Picked ", length(shuffle), " features: ", paste(shuffle, collapse = ", "))
  surv_cols <- colnames(full_data)[1:2]
  cmatrix <- full_data[, c(surv_cols, shuffle), drop = FALSE]
  return(cmatrix)
}

#' Regression call with timeout and retry logic
#'
#' @description Internal helper function used by \code{rebootRegression()}.
#' Executes regression on a bootstrap sample with a timeout limit.
#' If execution exceeds the time limit, a new subsample is drawn and the regression is retried once.
#'
#' @param cmatrix A data.frame with survival data and features.
#' @param group_size Integer. Number of features used in subsampling.
#' @param full_data Original full dataset used for resampling if needed.
#'
#' @return Named numeric vector of regression coefficients, or NULL if failure.
#' @keywords internal
reboot_regcall <- function(cmatrix, group_size, full_data) {
  
  # Double tryCatch to avoid calling the function itself
  result <- tryCatch(
    R.utils::withTimeout(reboot_regression(cmatrix), timeout = 60, onTimeout = "error"),
    error = function(e) {
      message("Regression timeout/error: retrying with new subsample...")
      new_matrix <- reboot_subsample(full_data, group_size)
      
      # Retries one more time
      tryCatch(
        R.utils::withTimeout(reboot_regression(new_matrix), timeout = 60, onTimeout = "error"),
        error = function(e2) {
          message("Retry failed: ", conditionMessage(e2))
          return(NULL)
        }
      )
    }
  )
  
  return(result)
}

#' Penalized (LASSO) Cox regression core for REBOOT pipeline
#'
#' @description Internal helper function used by \code{rebootRegression()}.
#' Fits a penalized Cox regression model using L1 (LASSO) penalty and returns coefficients.
#'
#' @param cmatrix A data.frame with survival data and features.
#'
#' @return Named numeric vector of regression coefficients.
#' @keywords internal
reboot_regression <- function(cmatrix) {
  
  # Ensures data has OS and OS.time data
  if (!all(c("OS", "OS.time") %in% colnames(cmatrix))) {stop("'cmatrix' must contain 'OS' and 'OS.time' names")}
  
  # Ensures survival formula exists correctly
  form <- survival::Surv(OS.time, OS) ~ .
  
  # Model exploration (initial cross-validation using fold = 10)
  fit1 <- tryCatch(
    penalized::profL1(form, data = cmatrix, fold = 10, plot = FALSE, trace = FALSE),
    error = function(e) {stop("'profL1' failed: ", conditionMessage(e))}
  )
  
  # Selection of the best lambda
  opt1 <- tryCatch(
    penalized::optL1(form, data = cmatrix, fold = fit1$fold, trace = FALSE),
    error = function(e) {stop("'optL1' failed: ", conditionMessage(e))}
  )
  
  # Final penalized model
  fit <- tryCatch(
    penalized::penalized(
      form,
      data = cmatrix,
      lambda1 = opt1$lambda,
      trace = FALSE
    ),
    error = function(e) {stop("'penalized' model failed: ", conditionMessage(e))}
  )
  
  coefs <- penalized::coefficients(fit, "all")
  if (is.null(coefs) || length(coefs) == 0) {stop("No coefficients returned by penalized regression model")}
  return(coefs)
}

###### MAIN (SURVIVAL - MODULE II) ######

#' Determine optimal ROC-based cutoff for REBOOT scores
#'
#' @description Internal helper function used by \code{rebootSurvival()} to determine the optimal cutoff for continuous
#' REBOOT scores using a time-dependent ROC curve and the Youden index criterion.
#'
#' The function automatically determines the ROC direction based on the AUC value estimated
#' with \code{survivalROC::survivalROC()} unless explicitly provided by the user.
#'
#' @param data A data.frame containing at least the columns \code{score}, \code{OS}, and \code{OS.time}.
#' @param direction Character. ROC direction passed to \code{OptimalCutpoints::optimal.cutpoints()}.
#' If \code{NULL} (default), direction is automatically inferred from the AUC value.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{cutoff}: numeric optimal cutoff value
#'   \item \code{auc}: numeric AUC value from the time-dependent ROC curve
#'   \item \code{optimal_cutpoint}: complete output from \code{OptimalCutpoints::optimal.cutpoints()}
#'   \item \code{direction}: ROC direction used
#' }
#'
#' @keywords internal
reboot_cutoff_roc <- function(data, direction = NULL)
{
  # Checks input type
  if (!is.data.frame(data)) {stop("Argument 'data' must be a data.frame")}
  
  # Validates input
  required_cols <- c("score", "OS", "OS.time")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {stop("Missing required columns: ", paste(missing_cols, collapse = ", "))}
  
  # Checks direction
  if (!is.null(direction) && (!is.character(direction) || length(direction) != 1 || !direction %in% c("<", ">")))
  {
    stop("Argument 'direction' must be NULL, '<', or '>'")
  }
  
  # Determines ROC direction automatically if not provided
  if (is.null(direction))
  {
    # Calculates ROC object
    auc_obj <- survivalROC::survivalROC(Stime = data$OS.time, status = data$OS, marker = data$score,
                                        predict.time = stats::median(data$OS.time, na.rm = TRUE),
                                        method = "NNE", span = 0.25*nrow(data)^(-0.20))
    
    # Extracts preliminary AUC used only for direction inference
    auc_preliminary <- as.numeric(auc_obj$AUC)
    
    # Gets direction
    direction <- ifelse(auc_preliminary < 0.5, ">", "<")
  }
  
  # Calculates optimal cutoff
  optimal_cutpoint <- OptimalCutpoints::optimal.cutpoints(X = "score", status = "OS", tag.healthy = 0,
                                                          methods = "Youden", data = data,
                                                          control = OptimalCutpoints::control.cutpoints(),
                                                          direction = direction)
  
  # Extracts final ROC AUC from optimal cutpoint analysis
  auc <- as.numeric(optimal_cutpoint$Youden$Global$measures.acc$AUC[["AUC"]])
  
  # Extracts cutoff
  cutoff <- as.numeric(optimal_cutpoint$Youden$Global$optimal.cutoff$cutoff)
  
  # Handles rare multiple-cutoff outputs
  if (length(cutoff) >= 2) {cutoff <- cutoff[1]}
  
  # Returns structured result
  return(list(cutoff = cutoff, auc = auc, optimal_cutpoint = optimal_cutpoint, direction = direction))
}

#' Run univariate Cox regression model for REBOOT scores
#'
#' @description Internal helper used by \code{rebootSurvival()} to perform univariate survival analysis using
#' categorized REBOOT scores.
#'
#' The function:
#' \itemize{
#'   \item Fits a univariate Cox proportional hazards model
#'   \item Tests proportional hazards assumptions
#'   \item Computes hazard ratios and confidence intervals
#'   \item Calculates median survival statistics
#'   \item Organizes all results into a structured list
#' }
#'
#' Input data must already contain:
#' \itemize{
#'   \item \code{OS}
#'   \item \code{OS.time}
#'   \item \code{score_group}
#' }
#'
#' @param data A data.frame containing survival information and categorized (binary) REBOOT scores.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{table}: formatted univariate Cox regression results table
#'   \item \code{model}: fitted \code{survival::coxph} object
#'   \item \code{fit}: fitted \code{survival::survfit} object
#'   \item \code{ph}: proportional hazards assumptions test results
#' }
#'
#' Returns \code{NULL} if model fitting fails.
#'
#' @keywords internal
reboot_uniCox_model <- function(data)
{
  # Checks input type
  if (!is.data.frame(data)) {stop("Argument 'data' must be a data.frame")}
  
  # Validates input
  required_cols <- c("OS", "OS.time", "score_group")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {stop("Missing required columns: ", paste(missing_cols, collapse = ", "))}
  
  # Runs analysis
  tryCatch({
    # Fits univariate Cox model
    uni_model <- survival::coxph(survival::Surv(OS.time, OS) ~ score_group, data = data)
    
    # Tests proportional hazards assumptions
    ph_test <- reboot_schoenfeld(model_object = uni_model, covariates = "score_group", is_multi = FALSE)
    
    # Computes hazard ratio statistics
    hr_model <- reboot_hazard_ratio(x = data$score_group, surv.time = data$OS.time, surv.event = data$OS,
                                    alpha = 0.05, method.test = "logrank", na.rm = TRUE)
    
    # Extracts statistics
    hazard_ratio <- paste0(round(hr_model$hazard.ratio, 4), " (95% CI, ", round(hr_model$lower, 4),
                           " - ", round(hr_model$upper, 4), ")")
    coef <- hr_model$coef
    logrank_pvalue <- hr_model$p.value
    
    # Counts samples per group
    low_high_samples <- paste0(nrow(data[data$score_group == "low", ]), "/", nrow(data[data$score_group == "high", ]))
    
    # Fits Kaplan-Meier estimator
    fit <- survival::survfit(survival::Surv(OS.time, OS) ~ score_group, data = data)
    
    # Extracts median survival statistics
    fit_summary <- summary(fit)$table
    median_survival_low <- paste0(fit_summary[2, 7], " (95% CI, ", fit_summary[2, 8], " - ", fit_summary[2, 9], ")")
    median_survival_high <- paste0(fit_summary[1, 7], " (95% CI, ", fit_summary[1, 8], " - ", fit_summary[1, 9], ")")
    
    # Assigns prognosis label
    if (!is.na(logrank_pvalue) && logrank_pvalue < 0.05 && coef < 0)
    {
      prognosis <- "better"
    } else if (!is.na(logrank_pvalue) && logrank_pvalue < 0.05 && coef > 0) {
      prognosis <- "worse"
    } else {prognosis <- "---"}
    
    # Builds results table
    result_table <- data.frame(feature = "score", control = "high", condition = "low", coefficient = round(coef, 4),
                               hazard.ratio = hazard_ratio, log.rank.pvalue = round(logrank_pvalue, 4),
                               low.high.samples = low_high_samples, median.survival.low = median_survival_low,
                               median.survival.high = median_survival_high, prognosis = prognosis,
                               stringsAsFactors = FALSE)
    
    # Returns structured result
    return(list(table = result_table, model = uni_model, fit = fit, ph = ph_test))
  }, error = function(e) {
    warning("Univariate Cox model failed: ", conditionMessage(e))
    return(NULL)
  })
}

#' Run multiple univariate Cox regressions for clinical covariates
#'
#' @description Fits independent univariate Cox proportional hazards models for each provided clinical covariate.
#'
#' @param data Data.frame containing survival and binary clinical information.
#' Must contain columns \code{OS} and \code{OS.time}.
#' @param covariates Character. Vector containing clinical covariate names.
#' @param p.cutoff Numeric. P-value threshold used to define prognosis groups (default: \code{0.2}).
#'
#' @return A data.frame containing univariate Cox regression statistics for all tested covariates.
#'
#' @keywords internal
reboot_uniCox_test <- function(data, covariates, p.cutoff = 0.2)
{
  # Validates input data
  required_cols <- c("OS", "OS.time")
  if (!all(required_cols %in% colnames(data))) {stop("Input data must contain columns: 'OS' and 'OS.time'")}
  
  # Validates covariates
  if (!is.character(covariates) || length(covariates) < 1) {stop("'covariates' must be a character vector")}
  missing_covariates <- covariates[!covariates %in% colnames(data)]
  if (length(missing_covariates) > 0) {stop(paste("Covariates not found in data:",
                                                  paste(missing_covariates, collapse = ", ")))}
  
  # Validates p-value cutoff
  if (!is.numeric(p.cutoff) || length(p.cutoff) != 1 || p.cutoff <= 0 || p.cutoff >= 1)
  {
    stop("'p.cutoff' must be a numeric value between 0 and 1")
  }
  
  # Runs multiple univariate Cox models for each variable
  results <- lapply(covariates, function(variable)
  {
    formula <- stats::as.formula(paste("survival::Surv(OS.time, OS) ~", variable))
    model <- tryCatch(survival::coxph(formula, data = data), error = function(e) NULL)
    
    # Skips failed models
    if (is.null(model)) {return(NULL)}
    
    # Extracts statistics
    coef_df <- as.data.frame(stats::coef(summary(model)))
    ci_df <- as.data.frame(summary(model)$conf.int)
    res <- cbind(
      variable = rownames(coef_df),
      coefficient = coef_df$coef,
      hazard.ratio.raw = ci_df$`exp(coef)`,
      Cox.pvalue = coef_df$`Pr(>|z|)`,
      lower.ci = ci_df$`lower .95`,
      upper.ci = ci_df$`upper .95`
    )
    
    # Organizes results as a data.frame
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    
    # Formats statistics
    res$coefficient <- as.numeric(res$coefficient)
    res$hazard.ratio.raw <- as.numeric(res$hazard.ratio.raw)
    res$Cox.pvalue <- as.numeric(res$Cox.pvalue)
    res$lower.ci <- as.numeric(res$lower.ci)
    res$upper.ci <- as.numeric(res$upper.ci)
    
    # Creates formatted HR column
    res$hazard.ratio <- paste0(round(res$hazard.ratio.raw, 4), " (95% CI, ",
                               round(res$lower.ci, 4), " - ", round(res$upper.ci, 4), ")")
    
    # Formats p-values
    res$Cox.pvalue <- round(res$Cox.pvalue, 4)
    
    # Defines prognosis
    res$prognosis <- "---"
    res$prognosis[!is.na(res$Cox.pvalue) & res$Cox.pvalue < p.cutoff & res$coefficient < 0] <- "better"
    res$prognosis[!is.na(res$Cox.pvalue) & res$Cox.pvalue < p.cutoff & res$coefficient > 0] <- "worse"
    
    # Returns final columns
    res[, c("variable", "hazard.ratio", "Cox.pvalue", "prognosis")]
  })
  
  # Combines all results into a list
  results <- Filter(Negate(is.null), results)
  
  # Checks if any final result is valid
  if (length(results) == 0)
  {
    warning("No valid univariate Cox models could be fitted")
    return(NULL)
  }
  
  # Transforms the list into a new consolidated data.frame
  univ_res <- do.call(rbind, results)
  rownames(univ_res) <- NULL
  
  return(univ_res)
}

#' Run multivariate Cox regression with optional bootstrap resampling
#'
#' @description Performs multivariate Cox regression using selected clinical covariates and REBOOT score groups.
#' Optionally applies bootstrap resampling to identify stable covariates prior to fitting the final model.
#'
#' @param data Data.frame containing survival, REBOOT score and clinical information.
#' Must contain columns \code{OS}, \code{OS.time}, and \code{score_group}.
#' @param univ_result Data.frame generated by \code{reboot_uniCox_test()}.
#' @param logrank_result Data.frame generated by \code{reboot_uniCox_model()}.
#' @param covariates Character. Vector containing covariates selected for multivariate analysis.
#' @param all_covariates Character. Vector containing all available clinical covariates.
#' @param use.bootstrap Logical. Whether bootstrap resampling should be used before final model fitting (default: \code{TRUE}).
#' @param bootstrap Numeric. Number of bootstrap iterations (default: \code{1}).
#' @param bootstrap.freq Numeric. Minimum frequency to retain variables after bootstrap resampling (default: \code{0.25}).
#' @param multi.p.cutoff Numeric. P-value used to define prognosis groups in multivariate models (default: \code{0.05}).
#'
#' @return A list containing:
#' \describe{
#'   \item{table}{Formatted multivariate Cox regression results.}
#'   \item{model}{Final fitted Cox proportional hazards model.}
#'   \item{covariates}{Covariates retained in the final model.}
#'   \item{bootstrap.table}{Bootstrap frequency table (if applicable).}
#'   \item{bootstrap.used}{Logical indicating whether bootstrap was used.}
#' }
#'
#' @keywords internal
reboot_multiCox_test <- function(data, univ_result, logrank_result, covariates, all_covariates,
                                 use.bootstrap = TRUE, bootstrap = 1, bootstrap.freq = 0.25, multi.p.cutoff = 0.05)
{
  # Validates input data
  required_cols <- c("OS", "OS.time", "score_group")
  if (!is.data.frame(data)) {stop("'data' must be a data.frame")}
  if (!all(required_cols %in% colnames(data))) {stop("Input data must contain columns: 'OS', 'OS.time', and 'score_group'")}
  
  # Validates univariate results for clinical variables
  if (!is.data.frame(univ_result)) {stop("'univ_result' must be a data.frame")}
  required_univ_cols <- c("variable", "hazard.ratio", "Cox.pvalue", "prognosis")
  if (!all(required_univ_cols %in% colnames(univ_result))) {stop("'univ_result' is missing required columns")}
  
  # Validates log-rank results for score
  if (!is.data.frame(logrank_result)) {stop("'logrank_result' must be a data.frame")}
  required_logrank_cols <- c("feature", "coefficient", "hazard.ratio", "log.rank.pvalue", "prognosis")
  if (!all(required_logrank_cols %in% colnames(logrank_result))) {stop("'logrank_result' is missing required columns")}
  
  # Validates covariates
  if (!is.character(covariates) || length(covariates) < 1) {stop("'covariates' must be a character vector")}
  if (!is.character(all_covariates) || length(all_covariates) < 1) {stop("'all_covariates' must be a character vector")}
  
  # Validates bootstrap parameters
  if (!is.logical(use.bootstrap) || length(use.bootstrap) != 1) {stop("'use.bootstrap' must be either TRUE or FALSE")}
  if (!is.numeric(bootstrap) || length(bootstrap) != 1 || bootstrap < 1) {stop("'bootstrap' must be a positive numeric value")}
  if (!is.numeric(bootstrap.freq) || length(bootstrap.freq) != 1 || bootstrap.freq <= 0 || bootstrap.freq > 1)
  {
    stop("'bootstrap.freq' must be a numeric value between 0 and 1")
  }
  if (!is.numeric(multi.p.cutoff) || length(multi.p.cutoff) != 1 || multi.p.cutoff <= 0 || multi.p.cutoff >= 1)
  {
    stop("'multi.p.cutoff' must be a numeric value between 0 and 1")
  }
  
  # Harmonizes score naming between univariate and multivariate analyses
  logrank_result$feature[logrank_result$feature == "score"] <- "score_group"
  
  # Ensures "score_group" is included
  covariates <- unique(c("score_group", covariates))
  
  # Internal helper to extract Cox statistics
  extract_cox_results <- function(model)
  {
    coef_df <- as.data.frame(stats::coef(summary(model)))
    ci_df <- as.data.frame(summary(model)$conf.int)
    
    res <- cbind(variable = rownames(coef_df), coefficient = coef_df$coef, hazard.ratio.raw = ci_df$`exp(coef)`,
                 Cox.pvalue = coef_df$`Pr(>|z|)`, lower.ci = ci_df$`lower .95`, upper.ci = ci_df$`upper .95`)
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    
    # Converts numeric columns
    numeric_cols <- c("coefficient", "hazard.ratio.raw", "Cox.pvalue", "lower.ci", "upper.ci")
    res[numeric_cols] <- lapply(res[numeric_cols], as.numeric)
    
    # Creates formatted hazard ratio column
    res$hazard.ratio <- paste0(round(res$hazard.ratio.raw, 4), " (95% CI, ",
                               round(res$lower.ci, 4), " - ", round(res$upper.ci, 4), ")")
    
    # Formats p-values
    res$Cox.pvalue <- round(res$Cox.pvalue, 4)
    
    # Defines prognosis
    res$prognosis <- NA
    res$prognosis[!is.na(res$Cox.pvalue) & res$Cox.pvalue < multi.p.cutoff & res$coefficient < 0] <- "better"
    res$prognosis[!is.na(res$Cox.pvalue) & res$Cox.pvalue < multi.p.cutoff & res$coefficient > 0] <- "worse"
    res[, c("variable", "hazard.ratio", "Cox.pvalue", "prognosis")]
  }
  
  # Builds multivariate formula
  multi_formula <- stats::as.formula(paste("survival::Surv(OS.time, OS) ~", paste(covariates, collapse = " + ")))
  
  # Initializes basic bootstrap objects
  bootstrap_table <- NULL
  bootstrap_used <- FALSE
  
  # Optional bootstrap resampling
  if (use.bootstrap)
  {
    # Runs bootstrap resampling
    message(paste0("Running bootstrap resampling with ", bootstrap, " iterations..."))
    resampling_results <- reboot_resampling(data = data, covariates = covariates, all_covariates = all_covariates,
                                            bootstrap = bootstrap, bootstrap.freq = bootstrap.freq,
                                            boot.p.cutoff = multi.p.cutoff)
    
    # Extracts main results from iterations
    covariates <- resampling_results$covariates
    bootstrap_table <- resampling_results$table
    bootstrap_used <- resampling_results$bootstrap.used
    
    # Rebuilds formula after variable selection
    multi_formula <- stats::as.formula(paste("survival::Surv(OS.time, OS) ~", paste(covariates, collapse = " + ")))
  }
  
  # Fits final multivariate Cox model
  model <- survival::coxph(formula = multi_formula, data = data)
  
  # Extracts multivariate statistics
  multi_res <- extract_cox_results(model)
  
  # Merges univariate (score) and multivariate (clinics) results
  final <- merge(univ_result, multi_res, by = "variable", all = TRUE)
  
  # Replaces score_group univariate values
  score_idx <- grepl("score_group", final$variable)
  final[score_idx, "hazard.ratio.x"] <- as.character(logrank_result$hazard.ratio)
  final[score_idx, "Cox.pvalue.x"] <- logrank_result$log.rank.pvalue
  final[score_idx, "prognosis.x"] <- as.character(logrank_result$prognosis)
  
  # Replaces missing values
  final[is.na(final)] <- "---"
  
  # Renames columns
  colnames(final) <- c("variable", "univariate.hazard.ratio", "univariate.Cox.pvalue", "univariate.prognosis",
                       "multivariate.hazard.ratio", "multivariate.Cox.pvalue", "multivariate.prognosis")
  
  # Extracts reference conditions
  tmp_variables <- c()
  tmp_conditions <- c()
  
  reference_covariates <- c("score_group", all_covariates)
  
  for (var in final$variable)
  {
    for (ref_var in reference_covariates)
    {
      if (grepl(ref_var, var))
      {
        tmp_variables <- append(tmp_variables, ref_var)
        ref_group <- substr(x = var, start = nchar(ref_var) + 1, stop = nchar(var))
        tmp_conditions <- append(tmp_conditions, ref_group)
      }
    }
  }
  
  final$variable <- tmp_variables
  final$condition <- tmp_conditions
  
  final <- final[, c("variable", "condition", "univariate.hazard.ratio", "univariate.Cox.pvalue", "univariate.prognosis",
                     "multivariate.hazard.ratio", "multivariate.Cox.pvalue", "multivariate.prognosis")]
  
  rownames(final) <- NULL
  
  return(list(table = final, model = model, covariates = covariates,
              bootstrap.table = bootstrap_table, bootstrap.used = bootstrap_used))
}

###### SUPPORT (SURVIVAL - MODULE II) ######

#' Test proportional hazards assumptions for Cox models
#'
#' @description Internal helper used by \code{rebootSurvival()} and \code{reboot_uniCox_model()} to evaluate
#' proportional hazards assumptions using Schoenfeld residuals from \code{survival::cox.zph()}.
#'
#' The function identifies which covariates satisfy the proportional hazards assumptions and
#' returns the approved variables for downstream analyses, as well as the global test result.
#'
#' @param model_object A fitted Cox proportional hazards model generated by \code{survival::coxph()}.
#' @param covariates Character vector containing covariate names included in the model.
#' @param is_multi Logical. Whether the model corresponds to a multivariate analysis (default: FALSE).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{variables}: covariates passing PH assumptions
#'   \item \code{dropped}: covariates failing PH assumptions
#'   \item \code{global_p}: global proportional hazards p-value
#'   \item \code{table}: full \code{cox.zph()} results table
#'   \item \code{object}: original \code{cox.zph()} object
#' }
#'
#' @keywords internal
reboot_schoenfeld <- function(model_object, covariates, is_multi = FALSE)
{
  # Validates model
  if (is.null(model_object))
  {
    if (is_multi)
    {
      warning("Multivariate Cox model could not be fitted due to convergence issues.")
    } else {warning("Univariate Cox model could not be fitted due to convergence issues.")}
    return(NULL)
  }
  
  # Validates model object
  if (!inherits(model_object, "coxph")) {stop("'model_object' must be a valid survival::coxph model")}
  
  # Validates covariates
  if (!is.character(covariates) || length(covariates) < 1) {stop("'covariates' must be a non-empty character vector")}
  
  # Validates multivariate flag
  if (!is.logical(is_multi) || length(is_multi) != 1) {stop("'is_multi' must be either TRUE or FALSE")}
  
  # Runs Schoenfeld residuals test
  test_ph <- survival::cox.zph(fit = model_object, global = TRUE)
  
  # Adjusts covariate names for multivariate models
  if (is_multi)
  {
    tmp_covariates <- character()
    covariates <- unique(c("score_group", covariates))
    new_covariates <- colnames(test_ph$y)
    
    for (var in new_covariates)
    {
      matched <- covariates[sapply(covariates, function(x) {grepl(paste0("^", x), var)})]
      if (length(matched) > 0) {tmp_covariates <- c(tmp_covariates, matched[[1]])}
    }
  } else {tmp_covariates <- covariates}
  
  # Validates recovered covariate names
  expected_n <- nrow(test_ph$table) - 1
  if (length(tmp_covariates) != expected_n) {stop("Mismatch between recovered covariate names and Schoenfeld test output")}
  
  # Updates labels
  rownames(test_ph$table) <- c(tmp_covariates, "GLOBAL")
  colnames(test_ph$y) <- tmp_covariates
  
  # Extracts PH statistics
  tmp_df <- as.data.frame(test_ph$table)
  global_p <- tmp_df[nrow(tmp_df), "p"]
  pvalues <- tmp_df[1:(nrow(tmp_df) - 1), "p", drop = FALSE]
  variables <- rownames(pvalues[pvalues$p > 0.05, , drop = FALSE])
  dropped <- rownames(pvalues[pvalues$p <= 0.05, , drop = FALSE])
  
  # Evaluates global PH assumptions
  if (global_p <= 0.05)
  {
    warning(paste0("Global proportional hazards assumptions not met (Global p = ", round(global_p, 4), ")."))
  } else {message(paste0("Global proportional hazards assumptions met (Global p = ", round(global_p, 4), ")."))}
  
  # Reports dropped covariates
  if (length(dropped) > 0)
  {
    warning(paste0("The following covariates failed proportional hazards assumptions and will be removed: ",
                   paste(dropped, collapse = ", ")))
  }
  
  # Returns structured output
  return(list(variables = variables, dropped = dropped, global_p = global_p, table = tmp_df, object = test_ph))
}

#' Extract hazard ratio statistics from a univariate Cox model
#'
#' @description Internal helper used by \code{reboot_uniCox_model()} to fit a univariate Cox proportional
#' hazards model and extract hazard ratio statistics, confidence intervals, coefficients, and p-values.
#'
#' @param x Predictor variable.
#' @param surv.time Numeric vector containing follow-up times.
#' @param surv.event Numeric or logical vector containing survival events.
#' @param alpha Numeric. Significance level used for confidence intervals (default: \code{0.05}).
#' @param method.test Character. Statistical test used to compute p-values.
#' One of \code{"logrank"}, \code{"likelihood.ratio"}, or \code{"wald"}.
#' @param na.rm Logical. Whether to remove incomplete observations (default: \code{FALSE}).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{hazard.ratio}: Hazard ratio
#'   \item \code{coef}: Cox regression coefficient
#'   \item \code{se}: Standard error
#'   \item \code{lower}: Lower confidence interval
#'   \item \code{upper}: Upper confidence interval
#'   \item \code{p.value}: Model p-value
#'   \item \code{n}: Number of observations
#'   \item \code{model}: Fitted Cox model
#' }
#'
#' Returns \code{NULL} if the Cox model cannot be fitted.
#'
#' @keywords internal
reboot_hazard_ratio <- function(x, surv.time, surv.event, alpha = 0.05,
                                method.test = c("logrank", "likelihood.ratio", "wald"), na.rm = FALSE)
{
  # Validates x
  if (length(x) != length(surv.time) || length(x) != length(surv.event))
  {
    stop("Arguments 'x', 'surv.time', and 'surv.event' must have the same length")
  }
  
  # Validates survival times
  if (!is.numeric(surv.time) || any(surv.time < 0, na.rm = TRUE))
  {
    stop("'surv.time' must be a numeric vector containing non-negative values")
  }
  
  # Validates survival event
  if (!is.numeric(surv.event)) {stop("Argument 'surv.event' must be numeric")}
  unique_events <- sort(unique(stats::na.omit(surv.event)))
  if (!setequal(unique_events, c(0, 1))) {stop("Argument 'surv.event' must contain only 0 and 1 values")}
  
  # Validates alpha
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
  {
    stop("Argument 'alpha' must be a numeric value between 0 and 1")
  }
  
  # Validates na.rm
  if (!is.logical(na.rm) || length(na.rm) != 1) {stop("Argument 'na.rm' must be TRUE or FALSE")}
  
  # Matches statistical test argument
  method.test <- match.arg(method.test)
  
  # Builds temporary dataset
  df <- data.frame(x = x, surv.time = surv.time, surv.event = surv.event)
  
  # Removes missing values if requested
  if (na.rm) {df <- stats::na.omit(df)}
  
  # Validates remaining observations
  if (nrow(df) < 3)
  {
    warning("Insufficient observations to fit Cox model")
    return(NULL)
  }
  
  # Fits Cox model safely
  model <- tryCatch(survival::coxph(survival::Surv(surv.time, surv.event) ~ x, data = df), error = function(e) NULL)
  
  # Stops if model fitting failed
  if (is.null(model)) {return(NULL)}
  
  # Extracts coefficients
  coef <- stats::coef(model)[1]
  se <- sqrt(diag(stats::vcov(model)))[1]
  
  # Selects statistical test
  stat <- switch(method.test, logrank = model$score, likelihood.ratio = 2 * diff(model$loglik), wald = model$wald.test)
  
  # Calculates p-value
  p.value <- stats::pchisq(stat, df = 1, lower.tail = FALSE)
  
  # Confidence interval multiplier
  z <- stats::qnorm(alpha / 2, lower.tail = FALSE)
  
  # Returns structured results
  return(list(hazard.ratio = exp(coef), coef = coef, se = se, lower = exp(coef - z * se),
              upper = exp(coef + z * se), p.value = p.value, n = model$n, model = model))
}

#' Perform bootstrap resampling for multivariate Cox regression
#'
#' @description Runs bootstrap-based variable stability analysis for multivariate Cox proportional hazards regression.
#' Used in \code{reboot_multiCox_test()}
#'
#' @param data Data.frame containing survival and clinical variables.
#' Must contain columns \code{OS}, \code{OS.time}, and \code{score_group}.
#' @param covariates Character. Vector containing selected covariates.
#' @param all_covariates Character. Vector containing all available clinical covariates.
#' @param bootstrap Numeric. Number of bootstrap iterations (default: \code{100}).
#' @param bootstrap.freq Numeric. Minimum frequency to retain variables after bootstrap resampling (default: \code{0.25}).
#' @param boot.p.cutoff Numeric. P-value threshold used to define prognosis groups (default: \code{0.05}).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{table}: Bootstrap frequency table.
#'   \item \code{covariates}: Selected covariates after bootstrap filtering.
#'   \item \code{bootstrap.used}: Logical indicating whether bootstrap was successfully used.
#' }
#'
#' @keywords internal
reboot_resampling <- function(data, covariates, all_covariates,
                              bootstrap = 100, bootstrap.freq = 0.25, boot.p.cutoff = 0.05)
{
  # Validates input data
  required_cols <- c("OS", "OS.time", "score_group")
  if (!all(required_cols %in% colnames(data))) {stop("Input data must contain columns: 'OS' and 'OS.time'")}
  
  # Validates covariates
  if (!is.character(covariates) || length(covariates) < 1) {stop("'covariates' must be a character vector")}
  if (!is.character(all_covariates) || length(all_covariates) < 1) {stop("'all_covariates' must be a character vector")}
  
  # Validates bootstrap iterations
  if (!is.numeric(bootstrap) || length(bootstrap) != 1 || bootstrap < 1) {stop("'bootstrap' must be a positive numeric value")}
  
  # Validates bootstrap frequency threshold
  if (!is.numeric(bootstrap.freq) || length(bootstrap.freq) != 1 || bootstrap.freq <= 0 || bootstrap.freq > 1)
  {
    stop("'bootstrap.freq' must be a numeric value between 0 and 1")
  }
  
  # Validates p-value cutoff
  if (!is.numeric(boot.p.cutoff) || length(boot.p.cutoff) != 1 || boot.p.cutoff <= 0 || boot.p.cutoff >= 1)
  {
    stop("'boot.p.cutoff' must be a numeric value between 0 and 1")
  }
  
  # Builds Cox formula
  multi_formula <- stats::as.formula(paste("survival::Surv(OS.time, OS) ~", paste(covariates, collapse = " + ")))
  
  # Internal helper to extract Cox results
  extract_boot_cox_results <- function(model, boot.p.cutoff)
  {
    coef_df <- as.data.frame(stats::coef(summary(model)))
    ci_df <- as.data.frame(summary(model)$conf.int)
    
    res <- cbind(variable = rownames(coef_df), coefficient = coef_df$coef, hazard.ratio.raw = ci_df$`exp(coef)`,
                 Cox.pvalue = coef_df$`Pr(>|z|)`, lower.ci = ci_df$`lower .95`, upper.ci = ci_df$`upper .95`)
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    
    res$coefficient <- as.numeric(res$coefficient)
    res$Cox.pvalue <- as.numeric(res$Cox.pvalue)
    res$prognosis <- NA
    res$prognosis[!is.na(res$Cox.pvalue) & res$Cox.pvalue < boot.p.cutoff & res$coefficient < 0] <- "better"
    res$prognosis[!is.na(res$Cox.pvalue) & res$Cox.pvalue < boot.p.cutoff & res$coefficient > 0] <- "worse"
    
    return(res[, c("variable", "prognosis")])
  }
  
  # Generates bootstrap samples
  bootstrap_data <- sjstats::bootstrap(data = data, n = bootstrap, size = 0.6)
  
  tables_list <- list()
  valid_iterations <- 0
  current_iteration <- 1
  
  while (valid_iterations < bootstrap && current_iteration <= length(bootstrap_data$strap))
  {
    boot_df <- reboot_bootstrap(raw_data = data, bootstrap_data = bootstrap_data, iteration = current_iteration, min.prop = 0.2)
    current_iteration <- current_iteration + 1
    
    # Skips invalid bootstrap samples
    if (is.null(boot_df)) {next}
    
    # Fits multivariate Cox model
    model <- tryCatch(survival::coxph(multi_formula, data = boot_df), error = function(e) NULL)
    if (is.null(model)) {next}
    
    # Extracts iteration results
    tables_list[[length(tables_list) + 1]] <- extract_boot_cox_results(model = model, boot.p.cutoff = boot.p.cutoff)
    valid_iterations <- valid_iterations + 1
  }
  
  # Ensures valid bootstrap results exist
  if (length(tables_list) == 0)
  {
    warning("Bootstrap resampling failed. Returning original covariates.")
    return(list(table = NULL, covariates = covariates, bootstrap.used = FALSE))
  }
  
  # Merges bootstrap results
  merged_table <- reboot_merge(tables_list = tables_list, covariates = all_covariates)
  
  # Selects frequent covariates
  frequency_cutoff <- bootstrap * bootstrap.freq
  selected_covariates <- merged_table$covariate[merged_table$frequency >= frequency_cutoff]
  
  # Ensures score_group remains included
  selected_covariates <- unique(c("score_group", selected_covariates))
  
  return(list(table = merged_table, covariates = selected_covariates, bootstrap.used = TRUE))
}

#' Validate bootstrap resampled datasets for multivariate Cox regression
#'
#' @description Validates whether covariates in a bootstrapped dataset remain binary and balanced for multivariate Cox regression.
#'
#' @param raw_data Data.frame used as source for bootstrap resampling.
#' @param bootstrap_data Bootstrap object generated by \code{sjstats::bootstrap()}.
#' @param iteration Numeric. Bootstrap iteration index.
#' @param min.prop Numeric. Minimum proportion for the least abundant category in binary variables (default: \code{0.2}).
#'
#' @return A bootstrap-resampled data.frame if valid, otherwise \code{NULL}.
#'
#' @keywords internal
reboot_bootstrap <- function(raw_data, bootstrap_data, iteration, min.prop = 0.2)
{
  # Validates input data
  if (!is.data.frame(raw_data)) {stop("'raw_data' must be a data.frame")}
  required_cols <- c("OS", "OS.time")
  if (!all(required_cols %in% colnames(raw_data))) {stop("'raw_data' must contain columns: 'OS' and 'OS.time'")}
  
  # Validates bootstrap data
  if (!is.list(bootstrap_data) || is.null(bootstrap_data$strap))
  {
    stop("'bootstrap_data' must be a valid bootstrap object generated by sjstats::bootstrap()")
  }
  
  # Validates iteration
  if (!is.numeric(iteration) || length(iteration) != 1 || iteration < 1) {stop("'iteration' must be a positive numeric value")}
  if (iteration > length(bootstrap_data$strap)) {stop("'iteration' exceeds the number of available bootstrap samples")}
  
  # Validates minimum proportion
  if (!is.numeric(min.prop) || length(min.prop) != 1 || min.prop <= 0 || min.prop >= 0.5)
  {
    stop("'min.prop' must be a numeric value between 0 and 0.5")
  }
  
  # Extracts bootstrap sample indices
  boot_indices <- bootstrap_data$strap[[iteration]]$id
  
  # Builds bootstrap dataset
  bootstrap_df <- as.data.frame(raw_data[boot_indices, , drop = FALSE])
  
  # Checks binary balance for covariates
  covariate_cols <- setdiff(colnames(bootstrap_df), c("OS", "OS.time"))
  for (col in covariate_cols)
  {
    current_var <- stats::na.omit(bootstrap_df[[col]])
    
    # Ensures binary variables
    if (length(current_var) == 0 || length(unique(current_var)) != 2) {return(NULL)}
    
    # Computes category proportions
    proportions <- prop.table(table(current_var))
    
    # Ensures minimum category proportion
    if (min(proportions) < min.prop) {return(NULL)}
  }
  
  return(bootstrap_df)
}

#' Merge bootstrap multivariate Cox regression results
#'
#' @description Summarizes bootstrap multivariate Cox regression results
#' by counting how frequently each covariate was selected across bootstrap iterations.
#'
#' @param tables_list List containing multivariate Cox regression result tables.
#' @param covariates Character vector containing original clinical covariates.
#'
#' @return A data.frame containing covariate selection frequencies across bootstrap iterations.
#'
#' @keywords internal
reboot_merge <- function(tables_list, covariates)
{
  # Validates tables list
  if (!is.list(tables_list) || length(tables_list) == 0) {stop("'tables_list' must be a non-empty list")}
  
  # Validates covariates
  if (!is.character(covariates) || length(covariates) < 1) {stop("'covariates' must be a character vector")}
  
  # Extracts significant variables from each iteration
  cox_tables <- lapply(tables_list, function(result_table)
  {
    if (!"prognosis" %in% colnames(result_table)) {stop("All tables in 'tables_list' must contain a 'prognosis' column")}
    result_table <- result_table[!is.na(result_table$prognosis), , drop = FALSE]
    if (nrow(result_table) == 0) {return(NULL)}
    dplyr::count(result_table, variable)
  })
  
  # Removes empty results
  cox_tables <- Filter(Negate(is.null), cox_tables)
  
  if (length(cox_tables) == 0)
  {
    warning("No valid bootstrap iterations available")
    return(NULL)
  }
  
  # Combines all bootstrap results
  tmp_table <- data.table::rbindlist(cox_tables, use.names = TRUE, fill = TRUE, idcol = "bootstrap.id")
  
  # Counts covariate frequencies
  final_table <- dplyr::count(tmp_table, variable)
  colnames(final_table) <- c("covariate", "frequency")
  
  # Matches covariate names to original variables
  original_covariates <- c("score_group", covariates)
  
  matched_covariates <- sapply(final_table$covariate, function(variable)
  {
    matched <- original_covariates[sapply(original_covariates, function(x) grepl(paste0("^", x), variable))]
    if (length(matched) > 0) matched[[1]] else variable
  })
  
  final_table$covariate <- matched_covariates
  
  return(final_table)
}
