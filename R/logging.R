# Internal environment for REBOOT
.reboot_env <- new.env(parent = emptyenv())

# Internal helper to generate timestamps
.timestamp <- function() {format(Sys.time(), "%Y-%m-%d %H:%M:%S")}

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
