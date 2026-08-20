# Internal environment for REBOOT
.reboot_env <- new.env(parent = emptyenv())

# Internal helper to generate timestamps
.timestamp <- function() {format(Sys.time(), "%Y-%m-%d %H:%M:%S")}

# Internal helper to validate a user-supplied 'seed' argument
.reboot_check_seed <- function(seed) {
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || seed %% 1 != 0)) {
    log_stop("Argument 'seed' must be a single integer or NULL")
  }
}

# Internal helper to snapshot the current global RNG state (if any) into .reboot_env,
# so it can be restored later via .reboot_restore_seed(). This allows REBOOT pipelines
# to set a local, reproducible seed without permanently altering the user's R session.
.reboot_save_seed <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    .reboot_env$old_seed <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    .reboot_env$old_seed <- NULL
  }
}

# Internal helper to restore the global RNG state previously saved with .reboot_save_seed().
# Intended to be registered with on.exit(..., add = TRUE) right after set.seed() is called,
# so the RNG stream is restored regardless of how the calling function exits (success or error).
.reboot_restore_seed <- function() {
  if (!is.null(.reboot_env$old_seed)) {
    assign(".Random.seed", .reboot_env$old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  .reboot_env$old_seed <- NULL
}

#' Log an informational message
#'
#' @description
#' Writes a message to the console and, if enabled, to the log file.
#'
#' @param ... Character inputs to be concatenated into a message.
#'
#' @keywords internal
log_message <- function(...) {
  msg <- paste(...)
  full_msg <- paste0("[", .timestamp(), "] [INFO] ", msg)
  message(full_msg)
  
  if (exists("log_con", envir = .reboot_env) && !is.null(.reboot_env$log_con) && isOpen(.reboot_env$log_con)) {
    writeLines(full_msg, .reboot_env$log_con)
  }
}

#' Log a warning message
#'
#' @description
#' Writes a warning to the console and, if enabled, to the log file.
#'
#' @param ... Character inputs to be concatenated into a warning message.
#'
#' @keywords internal
log_warning <- function(...) {
  msg <- paste(...)
  full_msg <- paste0("[", .timestamp(), "] [WARN] ", msg)
  warning(full_msg, call. = FALSE)
  
  if (exists("log_con", envir = .reboot_env) && !is.null(.reboot_env$log_con) && isOpen(.reboot_env$log_con)) {
    writeLines(full_msg, .reboot_env$log_con)
  }
}

#' Log an error and stop execution
#'
#' @description
#' Writes an error message to the log file (if enabled) and stops execution.
#'
#' @param ... Character inputs to be concatenated into an error message.
#'
#' @keywords internal
log_stop <- function(...) {
  msg <- paste(...)
  full_msg <- paste0("[", .timestamp(), "] [ERROR] ", msg)
  
  if (exists("log_con", envir = .reboot_env) && !is.null(.reboot_env$log_con) && isOpen(.reboot_env$log_con)) {
    writeLines(full_msg, .reboot_env$log_con)
  }
  
  stop(full_msg, call. = FALSE)
}
