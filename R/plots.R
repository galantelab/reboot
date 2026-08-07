utils::globalVariables(c("feature", "coefficient", "covariate", "frequency"))

#' Default ggplot2 theme for REBOOT visualizations
#'
#' @description Internal helper function that provides a consistent ggplot2 theme used across REBOOT plots
#' Can also be applied to custom plots for visual consistency.
#'
#' @return A \code{ggplot} theme object.
#' @export
reboot_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      text = ggplot2::element_text(face = "plain", colour = "black"),
      axis.text = ggplot2::element_text(face = "plain", colour = "black"),
      legend.text = ggplot2::element_text(face = "plain", colour = "black"),
      legend.title = ggplot2::element_text(face = "plain", colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black"),
      axis.line = ggplot2::element_line(colour = "black")
    )
}

###### REGRESSION (MODULE I) ######

#' Lollipop plot for REBOOT signature
#'
#' @description Generates a lollipop plot of top features ranked by coefficient.
#'
#' @param signature A data.frame with 'feature' and 'coefficient'.
#' @param outprefix Character. Output prefix for saving the plot (optional).
#' @param theme A ggplot2 theme object. Defaults to \code{reboot_theme()}.
#' @param width Numeric. Plot width in inches (Default: \code{8})
#' @param height Numeric. Plot height in inches (Default: \code{5})
#' @param dpi Numeric. Plot DPI resolution (Default: \code{300})
#'
#' @return A \code{ggplot} object.
#' @export
reboot_lollipop <- function(signature, outprefix = NULL, theme = reboot_theme(), width = 8, height = 5, dpi = 300) {
  
  message("Building lollipop plot...")
  
  # Validates input
  if (!all(c("feature", "coefficient") %in% colnames(signature))) {stop("'signature' must contain 'feature' and 'coefficient' columns")}
  
  # Validates theme
  if (!inherits(theme, "theme")) {stop("'theme' must be a valid ggplot2 theme object")}
  
  # Validates plot parameters
  if (!is.numeric(dpi) || length(dpi) != 1 || dpi <= 0) {stop("'dpi' must be a single positive numeric value")}
  if (!is.numeric(width) || length(width) != 1 || width <= 0) {stop("'width' must be a single positive numeric value")}
  if (!is.numeric(height) || length(height) != 1 || height <= 0) {stop("'height' must be a single positive numeric value")}
  
  # Prepares data
  tt <- signature[stats::complete.cases(signature), , drop = FALSE]
  tt <- tt[tt$coefficient != 0, , drop = FALSE]
  
  # Checks if signature is still valid
  if (nrow(tt) == 0) {stop("No valid coefficients available for plotting")}
  
  # Selects top 10 features
  tt <- dplyr::slice_max(tt, order_by = abs(coefficient), n = 10, with_ties = FALSE)
  
  # Orders features by coefficient values
  tt$feature <- forcats::fct_reorder(tt$feature, tt$coefficient)
  
  # Makes plot
  ll <- ggplot2::ggplot(tt, ggplot2::aes(x = feature, y = coefficient)) +
    ggplot2::geom_segment(ggplot2::aes(x = feature, xend = feature, y = 0, yend = coefficient),
                          linewidth = 1.2, color = "grey") +
    ggplot2::geom_point(size = 3, colour = "#1c9099") + theme +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  # Saves if requested
  if (!is.null(outprefix)) {
    if (!is.character(outprefix) || length(outprefix) != 1) {stop("'outprefix' must be a single character string")}
    fname <- paste0(outprefix, "_lollipop.pdf")
    ggplot2::ggsave(filename = fname, plot = ll, device = "pdf", width = width, height = height, dpi = dpi)
  }
  
  return(ll)
}

#' Histogram plot for REBOOT signature
#'
#' @description Generates a histogram of coefficient distribution from a REBOOT signature.
#'
#' @param signature A data.frame with 'feature' and 'coefficient'.
#' @param outprefix Character. Output prefix for saving the plot (optional).
#' @param theme A ggplot2 theme object. Defaults to \code{reboot_theme()}.
#' @param width Numeric. Plot width in inches (Default: \code{8})
#' @param height Numeric. Plot height in inches (Default: \code{5})
#' @param dpi Numeric. Plot DPI resolution (Default: \code{300})
#'
#' @return A \code{ggplot} object.
#' @export
reboot_histogram <- function(signature, outprefix = NULL, theme = reboot_theme(), width = 8, height = 5, dpi = 300) {
  
  message("Building histogram plot...")
  
  # Validates input
  if (!all(c("feature", "coefficient") %in% colnames(signature))) {stop("'signature' must contain 'feature' and 'coefficient' columns")}
  
  # Validates theme
  if (!inherits(theme, "theme")) {stop("'theme' must be a valid ggplot2 theme object")}
  
  # Validates plot sizes
  if (!is.numeric(dpi) || length(dpi) != 1 || dpi <= 0) {stop("'dpi' must be a single positive numeric value")}
  if (!is.numeric(width) || length(width) != 1 || width <= 0) {stop("'width' must be a single positive numeric value")}
  if (!is.numeric(height) || length(height) != 1 || height <= 0) {stop("'height' must be a single positive numeric value")}
  
  # Prepares data
  tt <- signature[stats::complete.cases(signature), , drop = FALSE]
  tt <- tt[tt$coefficient != 0, , drop = FALSE]
  
  # Checks if signature is still valid
  if (nrow(tt) == 0) {stop("No valid coefficients available for plotting")}
  
  # Makes plot
  pl <- ggplot2::ggplot(tt, ggplot2::aes(x = coefficient)) +
    ggplot2::geom_histogram(binwidth = 0.005, alpha = 1) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Coefficient", y = NULL) + theme
  
  # Saves if requested
  if (!is.null(outprefix)) {
    if (!is.character(outprefix) || length(outprefix) != 1) {stop("'outprefix' must be a single character string")}
    fname <- paste0(outprefix, "_histogram.pdf")
    ggplot2::ggsave(filename = fname, plot = pl, device = "pdf", width = width, height = height, dpi = dpi)
  }
  
  return(pl)
}

###### SURVIVAL (MODULE II) ######

#' Generate ROC curve plot for REBOOT scores
#'
#' @description Generates a time-dependent ROC curve plot
#'
#' @param x Either:
#' \itemize{
#'   \item The full list returned by \code{reboot_cutoff_roc()}
#'   \item An \code{optimal.cutpoints} object
#' }
#' @param outprefix Character. Output prefix for saving the plot (optional).
#' If \code{NULL}, the plot is not written to disk.
#' @param width Numeric. Plot width in inches (Default: \code{7}).
#' @param height Numeric. Plot height in inches (Default: \code{6}).
#'
#' @return A ROC curve generated by \code{OptimalCutpoints::plot.optimal.cutpoints()}.
#' @export
reboot_roc <- function(x, outprefix = NULL, width = 7, height = 6) {
  
  message("Building ROC curve...")
  
  # Extracts optimal cutpoint object
  if (is.list(x) && "optimal_cutpoint" %in% names(x)) {roc_obj <- x$optimal_cutpoint} else {roc_obj <- x}
  
  # Validates plot sizes
  if (!is.numeric(width) || length(width) != 1 || width <= 0) {stop("'width' must be a single positive numeric value")}
  if (!is.numeric(height) || length(height) != 1 || height <= 0) {stop("'height' must be a single positive numeric value")}
  
  # Opens PDF device if requested
  if (!is.null(outprefix))
  {
    if (!is.character(outprefix) || length(outprefix) != 1) {stop("'outprefix' must be a single character string")}
    grDevices::pdf(file = paste0(outprefix, "_roc_curve.pdf"), width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  # Creates ROC plot (which = 1), but not the PROC plot (which = 2)
  OptimalCutpoints::plot.optimal.cutpoints(x = roc_obj, legend = TRUE, which = 1, col = "blue", bg = "white")
  
  # Captures base R plot object
  plot_obj <- grDevices::recordPlot()
  
  return(plot_obj)
}

#' Generate Schoenfeld residual plots for REBOOT Cox models
#'
#' @description Creates Schoenfeld residual plots from proportional hazards assumption tests.
#' For multivariate Cox models, the custom internal \code{reboot_ggcoxzph()} is used to additionally display
#' global and individual Schoenfeld test p-values.
#'
#' @param x Either:
#' \itemize{
#'   \item A \code{cox.zph} object
#'   \item A list returned by \code{reboot_schoenfeld()}
#' }
#' @param is_multi Logical. Whether the model corresponds to a multivariate Cox analysis (default: \code{FALSE}).
#' @param outprefix Character. Output prefix for saving the plot (optional).
#' If \code{NULL}, the plot is not written to disk.
#' @param theme A ggplot2 theme object. Defaults to \code{reboot_theme()}.
#' @param width Numeric. Plot width in inches (optional). If \code{NULL}, the width is set automatically.
#' @param height Numeric. Plot height in inches (optional). If \code{NULL}, the height is set automatically.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return a list of \code{ggplot} objects.
#' @export
reboot_phplot <- function(x, is_multi = FALSE, outprefix = NULL,
                          theme = reboot_theme(), width = NULL, height = NULL, ...) {
  
  message("Building proportional hazards plot...")
  
  # Extracts cox.zph object
  if (inherits(x, "cox.zph"))
  {
    ph_object <- x
  } else if (is.list(x) && "object" %in% names(x)) {
    ph_object <- x$object
  } else {stop("Input must be a 'cox.zph' object or a valid reboot_schoenfeld() result")}
  
  # Validates theme
  if (!inherits(theme, "theme")) {stop("'theme' must be a valid ggplot2 theme object")}
  
  # Validates arguments
  if (!is.logical(is_multi) || length(is_multi) != 1) {stop("'is_multi' must be TRUE or FALSE")}
  
  # Validates plot sizes
  if (!is.null(width))
  {
    if (!is.numeric(width) || length(width) != 1 || width <= 0)
    {
      stop("'width' must be NULL or a single positive numeric value")
    }
  }
  
  if (!is.null(height))
  {
    if (!is.numeric(height) || length(height) != 1 || height <= 0)
    {
      stop("'height' must be NULL or a single positive numeric value")
    }
  }
  
  # Generates plots
  if (is_multi)
  {
    ph_plots <- reboot_ggcoxzph(fit = ph_object, ggtheme = theme, ...)
  } else {ph_plots <- survminer::ggcoxzph(fit = ph_object, ggtheme = theme, ...)}
  
  # Renames variables in plot labels
  ph_plots <- lapply(ph_plots, function(p)
  {
    if (!is.null(p$labels$y)) {p$labels$y <- gsub("score_group", "Score", p$labels$y)}
    return(p)
  })
  
  # Extracts number of plots
  n_plots <- length(ph_plots)
  
  # Dynamically defines number of columns
  if (n_plots <= 2) {ncol_plot <- n_plots}
  else if (n_plots <= 4) {ncol_plot <- 2}
  else if (n_plots <= 9) {ncol_plot <- 3}
  else {ncol_plot <- 4}
  
  # Dynamically defines number of rows
  nrow_plot <- ceiling(n_plots / ncol_plot)
  
  # Automatically defines figure dimensions if not provided
  if (is.null(width)) {width <- 4 * ncol_plot}
  if (is.null(height)) {height <- 4 * nrow_plot}
  
  # Combines plots into a single patchwork object
  ph_plots <- patchwork::wrap_plots(ph_plots, ncol = ncol_plot)
  
  # Saves plots if requested
  if (!is.null(outprefix))
  {
    if (!is.character(outprefix) || length(outprefix) != 1) {stop("'outprefix' must be a single character string")}
    filename <- paste0(outprefix, "_ph_assumptions_plot.pdf")
    grDevices::pdf(file = filename, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(ph_plots)
  }
  
  return(ph_plots)
}

#' Customized Schoenfeld residual plots for REBOOT multivariate models
#'
#' @description Internal plotting helper used by \code{reboot_phplot()} to generate custom Schoenfeld residual plots
#'
#' This function is based on \code{survminer::ggcoxzph()}, with modifications to:
#' \itemize{
#'   \item Include global proportional hazards test p-values in plot titles
#'   \item Improve visualization for multivariate Cox regression models
#'   \item Return plot lists compatible with REBOOT plotting workflows
#' }
#'
#' @param fit A \code{cox.zph} object generated by \code{survival::cox.zph()}.
#' @param resid Logical. Whether to display Schoenfeld residual points (default: \code{TRUE}).
#' @param se Logical. Whether to display confidence bands (default: \code{TRUE}).
#' @param df Numeric. Degrees of freedom used in spline smoothing (default: \code{4}).
#' @param nsmo Numeric. Number of spline smoothing points (default: \code{40}).
#' @param var Variables to plot. Can be numeric indices or variable names.
#' @param point.col Point color for residuals (default: \code{"red"}).
#' @param point.size Numeric. Point size (default: \code{1}).
#' @param point.shape Numeric. Point shape (default: \code{19}).
#' @param point.alpha Numeric. Point transparency (default: \code{1}).
#' @param caption Optional plot caption.
#' @param ggtheme ggplot2 theme object (default: \code{survminer::theme_survminer()}).
#' @param ... Additional arguments passed to \code{ggpubr::ggpar()}.
#'
#' @return A list of \code{ggplot} objects of class \code{"ggcoxzph"}.
#' @keywords internal
reboot_ggcoxzph <- function(fit, resid = T, se = T, df = 4, nsmo = 40, var, point.col = "red", point.size = 1,
                            point.shape = 19, point.alpha = 1, caption = NULL,
                            ggtheme = survminer::theme_survminer(), ...) {
  x <- fit
  if (!methods::is(x, "cox.zph")) {stop("Can't handle an object of class ", class(x))}
  
  xx <- x$x
  yy <- x$y
  d <- nrow(yy)
  df <- max(df)
  nvar <- ncol(yy)
  pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
  temp <- c(pred.x, xx)
  lmat <- splines::ns(temp, df = df, intercept = T)
  pmat <- lmat[1:nsmo, ]
  xmat <- lmat[-(1:nsmo), ]
  qmat <- qr(xmat)
  
  if (qmat$rank < df) {stop("Spline fit is singular, try a smaller degrees of freedom")}
  
  if (se)
  {
    bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
    xtx <- bk %*% t(bk)
    seval <- d * ((pmat %*% xtx) * pmat) %*% rep(1, df)
  }
  
  ylab <- paste("Beta(t) for", dimnames(yy)[[2]])
  
  if (missing(var))
  {
    var <- 1:nvar
  } else {
    if (is.character(var)) {var <- match(var, dimnames(yy)[[2]])}
    if (any(is.na(var)) || max(var) > nvar || min(var) < 1) {stop("Invalid variable requested")}
  }
  
  if (x$transform == "log")
  {
    xx <- exp(xx)
    pred.x <- exp(pred.x)
  } else if (x$transform != "identity") {
    xtime <- as.numeric(dimnames(yy)[[1]])
    indx <- !duplicated(xx)
    apr1 <- stats::approx(xx[indx], xtime[indx], seq(min(xx), max(xx), length = 17)[2 * (1:8)])
    temp <- signif(apr1$y, 2)
    apr2 <- stats::approx(xtime[indx], xx[indx], temp)
    xaxisval <- apr2$y
    xaxislab <- rep("", 8)
    for (i in 1:8) xaxislab[i] <- format(temp[i])
  }
  
  plots <- list()
  lapply(var, function(i) {
    invisible(round(x$table[i, 3],4) -> pval)
    invisible(round(x$table[nrow(x$table), 3],4) -> global)
    
    ggplot2::ggplot() + ggplot2::labs(title = paste0('Global Schoenfeld Test p: ', global),
                                      subtitle = paste0('Individual Schoenfeld Test p: ', pval)) +
      ggtheme + ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, vjust = .5, face = "bold",
                                                                  margin = ggplot2::margin(0, 0, 10, 0)),
                               plot.subtitle = ggplot2::element_text(hjust = 0, vjust = .5, face = "plain",
                                                                     margin = ggplot2::margin(10, 0, 10, 0))) -> gplot
    
    y <- yy[, i]
    yhat <- as.vector(pmat %*% qr.coef(qmat, y))
    
    if (resid) {yr <- range(yhat, y)} else {yr <- range(yhat)}
    
    if (se)
    {
      temp <- as.vector(2 * sqrt(x$var[i, i] * seval))
      yup <- yhat + temp
      ylow <- yhat - temp
      yr <- range(yr, yup, ylow)
    }
    
    if (x$transform == "identity") {
      gplot + ggplot2::geom_line(ggplot2::aes(x = pred.x, y = yhat)) + ggplot2::xlab("Time") +
        ggplot2::ylab(ylab[i]) + ggplot2::ylim(yr) -> gplot
    } else if (x$transform == "log") {
      gplot + ggplot2::geom_line(ggplot2::aes(x = log(pred.x), y = yhat)) + ggplot2::xlab("Time") +
        ggplot2::ylab(ylab[i]) + ggplot2::ylim(yr)  -> gplot
    } else {
      gplot + ggplot2::geom_line(ggplot2::aes(x = pred.x, y = yhat)) + ggplot2::xlab("Time") + ggplot2::ylab(ylab[i]) +
        ggplot2::scale_x_continuous(breaks = xaxisval, labels = xaxislab) + ggplot2::ylim(yr) -> gplot
    }
    
    if (resid)
    {
      gplot <- gplot + ggplot2::geom_point(ggplot2::aes(x = xx, y = y), col = point.col,
                                           shape = point.shape, size = point.size, alpha = point.alpha)
    }
    
    if (se)
    {
      gplot <- gplot + ggplot2::geom_line(ggplot2::aes(x = pred.x, y = yup), lty = "dashed") +
        ggplot2::geom_line(ggplot2::aes( x = pred.x, y = ylow), lty = "dashed")
    }
    
    ggpubr::ggpar(gplot, ...)
  }) -> plots
  
  names(plots) <- var
  class(plots) <- c("ggcoxzph", "ggsurv", "list")
  
  # case of multivariate Cox
  if ("GLOBAL" %in% rownames(x$table)) {global_p <- x$table["GLOBAL", 3]} else {global_p <- NULL}
  
  attr(plots, "global_pval") <- global_p
  attr(plots, "caption") <- caption
  plots
}

#' Generate Kaplan-Meier survival plot for REBOOT scores
#'
#' @description Creates a Kaplan-Meier survival curve using categorized REBOOT scores generated by
#' \code{rebootSurvival()}. Input data must already contain a binary \code{score_group} variable.
#'
#' @param data Data.frame containing survival information and categorized scores.
#' Must contain columns \code{OS}, \code{OS.time}, and \code{score_group}.
#' @param cutoff Numeric. Cutoff value used to categorize scores. Displayed in legend labels (optional).
#' @param outprefix Character. Output prefix for saving the plot (optional).
#' If \code{NULL}, the plot is not written to disk.
#' @param pval Logical, Numeric or Character. Value to display in the plot. Typically extracted from \code{reboot_uniCox_model()} (optional).
#' @param palette Character vector. Color palette for Kaplan-Meier curves (default: \code{c("#D73027", "#1A9850")}).
#' @param ylab Character. Title for the Y-axis (default: \code{"Survival probability"}).
#' @param xlab Character. Title for the X-axis (default: \code{"Time"}).
#' @param title Character. Title for the plot (default: \code{NULL}).
#' @param legend.title Character. Title for the legend (default: \code{"Reboot groups"}).
#' @param theme ggplot2 theme used in the Kaplan-Meier plot (default: \code{survminer::theme_survminer()}).
#' @param tables.theme ggplot2 theme used in the risk table (default: \code{survminer::theme_cleantable()}).
#' @param width Numeric. Plot width in inches (default: \code{8}).
#' @param height Numeric. Plot height in inches (default: \code{5}).
#'
#' @return A Kaplan-Meier plot generated by \code{survminer::ggsurvplot()}.
#' @export
reboot_km <- function(data, cutoff = NULL, outprefix = NULL, pval = FALSE, palette = c("#D73027", "#1A9850"),
                      ylab = "Survival probability", xlab = "Time", title = NULL, legend.title = "Reboot groups",
                      theme = survminer::theme_survminer(), tables.theme = survminer::theme_cleantable(),
                      width = 8, height = 5) {
  
  message("Building Kaplan-Meier curves...")
  
  # Validates input
  required_cols <- c("OS", "OS.time", "score_group")
  if (!all(required_cols %in% colnames(data))) {stop("Input data must contain columns: 'OS', 'OS.time', and 'score_group'")}
  
  # Validates p-value parameter
  if (!(is.logical(pval) || is.numeric(pval) || is.character(pval))) {stop("'pval' must be logical, numeric, or character")}
  if (is.logical(pval) && length(pval) != 1) {stop("'pval' must be a single logical value")}
  if (is.character(pval) && length(pval) != 1) {stop("'pval' must be a single character string")}
  if (is.numeric(pval))
  {
    if (length(pval) != 1 || pval < 0 || pval > 1) {stop("'pval' must be a single numeric value between 0 and 1")}
  }
  
  # Validates plot parameters
  if (!is.null(cutoff)) {if (!is.numeric(cutoff) || length(cutoff) != 1) {stop("'cutoff' must be a single numeric value")}}
  if (!is.character(palette) || length(palette) != 2) {stop("'palette' must be a character vector of length 2")}
  if (!is.character(ylab) || length(ylab) != 1) {stop("'ylab' must be a single character string")}
  if (!is.character(xlab) || length(xlab) != 1) {stop("'xlab' must be a single character string")}
  if (!is.null(title) && (!is.character(title) || length(title) != 1))
  {
    stop("'title' must be NULL or a single character string")
  }
  if (!is.null(legend.title) && (!is.character(legend.title) || length(legend.title) != 1))
  {
    stop("'legend.title' must be NULL or a single character string")
  }
  
  # Validates plot sizes
  if (!is.numeric(width) || length(width) != 1 || width <= 0) {stop("'width' must be a single positive numeric value")}
  if (!is.numeric(height) || length(height) != 1 || height <= 0) {stop("'height' must be a single positive numeric value")}
  
  # Ensures at least two groups exist
  if (length(unique(stats::na.omit(data$score_group))) < 2)
  {
    warning("Data could not be partitioned into low/high scores")
    return(NULL)
  }
  
  # Fits Kaplan-Meier model
  fit <- survival::survfit(survival::Surv(OS.time, OS) ~ score_group, data = data)
  
  # Creates legend labels
  if (!is.null(cutoff))
  {
    legend_labels <- c(paste0("score > ", round(cutoff, 4)), paste0("score <= ", round(cutoff, 4)))
  } else {legend_labels <- unique(as.character(stats::na.omit(data$score_group)))}
  
  # Generates Kaplan-Meier plot
  km_plot <- suppressMessages(suppressWarnings(survminer::ggsurvplot(
    fit = fit,
    data = data,
    risk.table = "abs_pct",
    surv.scale = "percent",
    palette = palette,
    ylab = ylab,
    ylim = c(0, 1),
    break.y.by = 0.25,
    tables.height = 0.2,
    risk.table.title = "Number at risk (%)",
    tables.y.text = FALSE,
    xlab = xlab,
    tables.theme = tables.theme,
    ggtheme = theme,
    conf.int = FALSE,
    linetype = 1,
    censor.shape = 73,
    censor = TRUE,
    censor.size = 4,
    pval = pval,
    font.legend = 14,
    font.x = 18,
    font.y = 18,
    font.tickslab = 14,
    pval.size = 5,
    pval.coord = c(0, 0.05),
    title = title,
    legend = c(0.8, 0.9),
    legend.title = legend.title,
    linewidth = 1,
    fontsize = 3,
    axes.offset = TRUE,
    surv.median.line = "none",
    legend.labs = legend_labels
  )))
  
  # Adjusts risk table formatting
  km_plot$table <- km_plot$table + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0),
                                                  plot.margin = ggplot2::margin(5, 5, 5, 30))
  
  # Saves plot
  if (!is.null(outprefix))
  {
    if (!is.character(outprefix) || length(outprefix) != 1) {stop("'outprefix' must be a single character string")}
    grDevices::pdf(file = paste0(outprefix, "_km_plot.pdf"), width = width, height = height, onefile = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    suppressMessages(suppressWarnings(print(km_plot)))
  }
  
  return(km_plot)
}

#' Generate forest plot for multivariate Cox regression model
#'
#' @description Creates a forest plot from a multivariate Cox proportional hazards model.
#'
#' @param model_object A fitted \code{survival::coxph} model object.
#' @param outprefix Character. Output prefix for saving the plot (optional).
#' @param panels List. Custom forestmodel panels. If \code{NULL}, default REBOOT panels are used (optional).
#' @param hr_breaks Numeric. Vector defining hazard ratio axis breaks (default: \code{c(0.25, 0.5, 1, 2, 4, 8)}).
#' @param width Numeric. Plot width in inches (default: \code{8}).
#' @param height Numeric. Plot height in inches (default: \code{5}).
#'
#' @return A forest plot object generated by \code{forestmodel::forest_model()}.
#' @export
reboot_forest <- function(model_object, outprefix = NULL, panels = NULL,
                          hr_breaks = c(0.25, 0.5, 1, 2, 4, 8), width = 8, height = 5) {
  
  message("Building forest plot...")
  
  # Validates model object
  if (!inherits(model_object, "coxph")) {stop("'model_object' must be a valid 'survival::coxph' object")}
  
  # Validates custom panel
  if (!is.null(panels) && !is.list(panels)) {stop("'panels' must be NULL or a list")}
  
  # Validates HR breaks
  if (!is.numeric(hr_breaks) || length(hr_breaks) < 1 || any(hr_breaks <= 0)) {stop("'hr_breaks' must contain positive numeric values")}
  
  # Validates plot sizes
  if (!is.numeric(width) || length(width) != 1 || width <= 0) {stop("'width' must be a single positive numeric value")}
  if (!is.numeric(height) || length(height) != 1 || height <= 0) {stop("'height' must be a single positive numeric value")}
  
  # Uses default custom forestmodel panels if none provided
  if (is.null(panels))
  {
    panels <- list(list(width = 0.01),
                   list(width = 0.1, display = ~ifelse(variable == "score_group", "Score", variable),
                        fontface = "plain", heading = "Variable"),
                   list(width = 0.1, display = ~level, fontface = "italic", heading = "Group"),
                   list(width = 0.05, display = ~n, hjust = 1, fontface = "plain", heading = "N"),
                   list(width = 0.01, item = "vline", hjust = 0.5),
                   list(width = 0.7, item = "forest", hjust = 0.5, heading = "Hazard Ratio", linetype = "dashed", line_x = 0),
                   list(width = 0.01, item = "vline", hjust = 0.5),
                   list(width = 0.05, fontface = "plain", heading = "HR (95% CI)",
                        display = ~ifelse(reference, "Reference",
                                          sprintf("%0.2f (%0.2f - %0.2f)",
                                                  trans(estimate), trans(conf.low), trans(conf.high))), display_na = NA),
                   list(width = 0.05, fontface = "plain",
                        display = ~ifelse(reference, "", ifelse(p.value < 0.0001, "<0.0001",
                                                                round(x = p.value, digits = 4))),
                        display_na = NA, hjust = 1, heading = "p-value"),
                   list(width = 0.01))
  }
  
  # Extracts model confidence intervals
  sum_model <- summary(model_object)
  conf_int <- sum_model$conf.int
  hr_min <- min(conf_int[, "lower .95"], na.rm = TRUE)
  hr_max <- max(conf_int[, "upper .95"], na.rm = TRUE)
  
  # Avoids invalid values
  hr_min <- max(hr_min, 1e-6)
  
  # Clips extreme values
  hr_min <- max(hr_min, 0.001)
  hr_max <- min(hr_max, 100)
  
  # Creates symmetric log scale around HR = 1
  max_range <- max(abs(log(hr_min)), abs(log(hr_max)))
  
  # Defines symmetrical boundaries
  max_range <- round(max_range, 2)
  limits_log <- c(-max_range, max_range)
  
  # Adjusts HR breaks to plotting range
  hr_breaks <- hr_breaks[hr_breaks >= exp(limits_log[1]) & hr_breaks <= exp(limits_log[2])]
  if (length(hr_breaks) < 2) {hr_breaks <- c(0.5, 1, 2)}
  breaks_log <- log(hr_breaks)
  
  # Saves plot
  if (!is.null(outprefix))
  {
    if (!is.character(outprefix) || length(outprefix) != 1) {stop("'outprefix' must be a single character string")}
    grDevices::pdf(file = paste0(outprefix, "_forest_plot.pdf"), width = width, height = height, onefile = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  # Generates forest plot
  forest_plot <- suppressMessages(suppressWarnings(forestmodel::forest_model(model = model_object, exponentiate = TRUE,
                                                                             panels = panels, breaks = breaks_log,
                                                                             limits = limits_log,
                                                                             factor_separate_line = FALSE,
                                                                             recalculate_width = TRUE,
                                                                             recalculate_height = TRUE)))
  if (!is.null(outprefix)) {suppressMessages(suppressWarnings(print(forest_plot)))}
  
  return(forest_plot)
}

#' Generate bootstrap covariate frequency barplot
#'
#' @description Creates a barplot with the frequency of covariates selected in multivariate Cox regression with bootstrap.
#'
#' @param plot_df Data.frame generated by \code{reboot_merge()} or extracted from \code{reboot_multiCox_test()}
#' bootstrap results. Must contain columns \code{covariate} and \code{frequency}.
#' @param outprefix Character. Output prefix for saving the plot (optional).
#' @param frequency.threshold Numeric. Frequency threshold displayed as a horizontal dashed line (default: \code{25}).
#' @param theme A ggplot2 theme object. Defaults to \code{reboot_theme()}.
#' @param bar.width Numeric. Width of bars (default: \code{0.5}).
#' @param ylab Character. Title for the Y-axis (default: \code{Frequency (\%)}).
#' @param xlab Character. Title for the X-axis (default: empty string).
#' @param width Numeric. Plot width in inches (default: \code{8}).
#' @param height Numeric. Plot height in inches (default: \code{5}).
#'
#' @return A \code{ggplot} object.
#' @export
reboot_barplot <- function(plot_df, outprefix = NULL, frequency.threshold = 25, theme = reboot_theme(),
                           bar.width = 0.5, ylab = "Frequency (%)", xlab = "", width = 8, height = 5)
{
  
  message("Building bar plot...")
  
  # Validates input data
  if (!is.data.frame(plot_df)) {stop("'plot_df' must be a data.frame")}
  required_cols <- c("covariate", "frequency")
  if (!all(required_cols %in% colnames(plot_df))) {stop("'plot_df' must contain columns: 'covariate' and 'frequency'")}
  
  # Validates frequency threshold
  if (!is.numeric(frequency.threshold) || length(frequency.threshold) != 1 ||
      frequency.threshold < 0 || frequency.threshold > 100)
  {
    stop("'frequency.threshold' must be a numeric value between 0 and 100")
  }
  
  # Validates bar width
  if (!is.numeric(bar.width) || length(bar.width) != 1 || bar.width <= 0) {stop("'bar.width' must be a single positive numeric value")}
  
  # Validates theme
  if (!inherits(theme, "theme")) {stop("'theme' must be a valid ggplot2 theme object")}
  
  # Validates labels for the axis
  if (!is.character(ylab) || length(ylab) != 1) {stop("'ylab' must be a single character string")}
  if (!is.character(xlab) || length(xlab) != 1) {stop("'xlab' must be a single character string")}
  
  # Validates plot sizes
  if (!is.numeric(width) || length(width) != 1 || width <= 0) {stop("'width' must be a single positive numeric value")}
  if (!is.numeric(height) || length(height) != 1 || height <= 0) {stop("'height' must be a single positive numeric value")}
  
  # Ensures frequency is numeric
  plot_df$frequency <- as.numeric(plot_df$frequency)
  if (any(is.na(plot_df$frequency))) {stop("'frequency' column contains NA values")}
  if (any(plot_df$frequency < 0 || plot_df$frequency > 100)) {stop("'frequency' values must be between 0 and 100")}
  
  # Sorts variables by decreasing frequency
  plot_df <- plot_df[order(plot_df$frequency, decreasing = TRUE), , drop = FALSE]
  
  # Generates barplot
  final_plot <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = stats::reorder(covariate, -frequency), y = frequency)) +
    ggplot2::geom_col(width = bar.width, color = "grey", fill = "grey") +
    ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
    ggplot2::geom_hline(yintercept = frequency.threshold, color = "#BB1136", linetype = "dashed", linewidth = 0.5) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
    theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  
  # Saves plot
  if (!is.null(outprefix))
  {
    if (!is.character(outprefix) || length(outprefix) != 1) {stop("'outprefix' must be a single character string")}
    grDevices::pdf(file = paste0(outprefix, "_bootstrap_barplot.pdf"), width = width, height = height, onefile = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(final_plot)
  }
  
  return(final_plot)
}
