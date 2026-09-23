#!/usr/bin/env Rscript

library(argparse)
library(Reboot)

parser <- ArgumentParser(
  prog = "Reboot.R",
  description = "Regression and survival tool with a multivariate bootstrap approach"
)

subparsers <- parser$add_subparsers(
  dest = "sub_name",
  metavar = "<subcommand>",
  help = "choose only one option"
)

parser$add_argument(
  '-v', '--version',
  action = 'version',
  version = 'reboot 1.2.0'
)

########## create the parser for the "REGRESSION" command ##########
parser_reg <- subparsers$add_parser(
  'regression',
  help = "generates a molecular signature through penalized LASSO Cox regression"
)

parser_reg$add_argument(
  "-I", "--filein",
  type = "character",
  dest = "fname",
  metavar = '',
  help = "Input file name. Tab separated values (tsv) file containing expression and survival parameters"
)

parser_reg$add_argument(
  "-O", "--outprefix",
  type = "character",
  dest = "out",
  metavar = '',
  default = "reboot",
  help = 'Output file prefix. Default: reboot'
)

parser_reg$add_argument(
  "-B","--bootstrap",
  type = "integer",
  dest = "boot_iter",
  metavar = '',
  default = 1,
  help = "Number of iterations for bootstrap simulation (int). Default: 1"
)

parser_reg$add_argument(
  "-G", "--groupsize",
  type = "integer",
  dest = "nel",
  metavar = '',
  default = 10,
  help = "Number of genes/transcripts for each bootstrap simulation (int). Default: 10"
)

parser_reg$add_argument(
  "-P", "--pcentfilter",
  type = "double",
  dest = "pf",
  metavar = '',
  default = 0.3,
  help = "Percentage of correlated feature pairs allowed in each iteration (double). Default: 0.3"
)

parser_reg$add_argument(
  "-T", "--type",
  type = "character",
  dest = "type",
  metavar = '',
  choices = c("gene", "transcript"),
  default = "gene",
  help = 'Type of transcriptome data: gene or transcript (character). Default: gene'
)

parser_reg$add_argument(
  "-V", "--varfilter",
  type = "double",
  dest = "var",
  metavar = '',
  default = 0.01,
  help = 'Minimum normalized variance (0-1) required for each feature among samples (double). Default: 0.01'
)

parser_reg$add_argument(
  "-F", "--force",
  dest = "force",
  action = "store_true",
  default = FALSE,
  help = "Choose -F to bypass OS and OS.time filters"
)

parser_reg$add_argument(
  "--seed",
  type = "integer",
  dest = "seed",
  metavar = '',
  default = 123,
  help = "Random seed used for reproducible bootstrap regression (integer). Default: 123"
)

parser_reg$add_argument(
  "--ncores",
  type = "integer",
  dest = "ncores",
  metavar = '',
  default = 1,
  help = "Number of CPU cores used for parallel bootstrap regression (integer). Default: 1"
)

parser_reg$add_argument(
  "--followup",
  type = "double",
  dest = "followup",
  metavar = '',
  default = NULL,
  help = "Maximum followup time (double). Default: NULL"
)
#####################################################################################################################

########## create the parser for the "SURVIVAL" command ##########
parser_sur <- subparsers$add_parser(
  'survival',
  help = "applies a molecular signature in survival analysis for prognosis"
)

parser_sur$add_argument(
  "-I", "--filein",
  type = "character",
  dest = "fname",
  metavar = '',
  help = "Input file name. Tab separated values (tsv) file containing expression and survival parameters"
)

parser_sur$add_argument(
  "-O", "--outprefix",
  type = "character",
  dest = "out",
  metavar = '',
  default = "reboot",
  help = 'Output file prefix. Default: reboot'
)

parser_sur$add_argument(
  "-S", "--signature",
  type = "character",
  dest = "sig",
  metavar = '',
  help = 'Tab separated values (tsv) file containing a set of features and corresponding regression coefficients'
)

parser_sur$add_argument(
  "-M", "--multivariate",
  dest = "multi",
  action = "store_true",
  default = FALSE,
  help = 'If clinical variables should be included, choose -M. This option is tied with -C option'
)

parser_sur$add_argument(
  "-C", "--clinical",
  type = "character",
  dest = "clin_file",
  metavar = '',
  default = "",
  help = 'Tab separated values (tsv) file containing binary categorical variables only. Required if -M option is chosen'
)

parser_sur$add_argument(
  "-R", "--roc",
  dest = "roc_curve",
  action = "store_true",
  default = FALSE,
  help = 'If molecular score should be categorized according to a ROC curve instead of the median, choose -R'
)

parser_sur$add_argument(
  "-B","--bootstrap",
  type = "integer",
  dest = "boot_iter",
  metavar = '',
  default = 1,
  help = "Number of iterations for bootstrap simulation (int). Default: 1"
)

parser_sur$add_argument(
  "-V", "--varfilter",
  type = "double",
  dest = "var",
  metavar = '',
  default = 0.01,
  help = 'Minimum normalized variance (0-1) required for follow up time (double). Default: 0.01'
)

parser_sur$add_argument(
  "-F", "--force",
  dest = "force",
  action = "store_true",
  default = FALSE,
  help = "Choose -F to bypass OS and OS.time filters"
)

parser_sur$add_argument(
  "--p-cutoff",
  type = "double",
  dest = "p_cutoff",
  metavar = '',
  default = 0.2,
  help = "P-value threshold used to select clinical covariates from univariate Cox regression (double). Default: 0.2"
)

parser_sur$add_argument(
  "--followup",
  type = "double",
  dest = "followup",
  metavar = '',
  default = NULL,
  help = "Maximum followup time (double). Default: NULL"
)
#####################################################################################################################

########## create the parser for the "COMPLETE" command ##########
parser_all <- subparsers$add_parser(
  'complete',
  help = "generates and applies a molecular signature (integrated analysis)"
)

parser_all$add_argument(
  "-I", "--filein",
  type = "character",
  dest = "fname",
  metavar = '',
  help = "Input file name. Tab separated values (tsv) file containing expression and survival parameters"
)

parser_all$add_argument(
  "-O", "--outprefix",
  type = "character",
  dest = "out",
  metavar = '',
  default = "reboot",
  help = 'Output file prefix. Default: reboot'
)

parser_all$add_argument(
  "-M", "--multivariate",
  dest = "multi",
  action = "store_true",
  default = FALSE,
  help = 'If clinical variables should be included, choose -M. This option is tied with -C option'
)

parser_all$add_argument(
  "-C", "--clinical",
  type = "character",
  dest = "clin_file",
  metavar = '',
  default = "",
  help = 'Tab separated values (tsv) file containing binary categorical variables only. Required if -M option is chosen'
)

parser_all$add_argument(
  "-R", "--roc",
  dest = "roc_curve",
  action = "store_true",
  default = FALSE,
  help = 'If molecular score should be categorized according to a ROC curve instead of the median, choose -R'
)

parser_all$add_argument(
  "-B","--bootstrap",
  type = "integer",
  dest = "boot_iter",
  metavar = '',
  default = 1,
  help = "Number of iterations for bootstrap simulation (int). Default: 1"
)

parser_all$add_argument(
  "-G", "--groupsize",
  type = "integer",
  dest = "nel",
  metavar = '',
  default = 10,
  help = "Number of genes/transcripts for each bootstrap simulation (int). Default: 10"
)

parser_all$add_argument(
  "-P", "--pcentfilter",
  type = "double",
  dest = "pf",
  metavar = '',
  default = 0.3,
  help = "Percentage of correlated feature pairs allowed in each iteration (double). Default: 0.3"
)

parser_all$add_argument(
  "-T", "--type",
  type = "character",
  dest = "type",
  metavar = '',
  choices = c("gene", "transcript"),
  default = "gene",
  help = 'Type of transcriptome data: gene or transcript (character). Default: gene'
)

parser_all$add_argument(
  "-V", "--varfilter",
  type = "double",
  dest = "var",
  metavar = '',
  default = 0.01,
  help = 'Minimum normalized variance (0-1) required for each feature and follow up time (double). Default: 0.01'
)

parser_all$add_argument(
  "-F", "--force",
  dest = "force",
  action = "store_true",
  default = FALSE,
  help = "Choose -F to bypass OS and OS.time filters"
)

parser_all$add_argument(
  "--seed",
  type = "integer",
  dest = "seed",
  metavar = '',
  default = 123,
  help = "Random seed used for reproducible bootstrap regression (integer). Default: 123"
)

parser_all$add_argument(
  "--ncores",
  type = "integer",
  dest = "ncores",
  metavar = '',
  default = 1,
  help = "Number of CPU cores used for parallel bootstrap regression (integer). Default: 1"
)

parser_all$add_argument(
  "--followup",
  type = "double",
  dest = "followup",
  metavar = '',
  default = NULL,
  help = "Maximum followup time (double). Default: NULL"
)

parser_all$add_argument(
  "--p-cutoff",
  type = "double",
  dest = "p_cutoff",
  metavar = '',
  default = 0.2,
  help = "P-value threshold used to select clinical covariates from univariate Cox regression (double). Default: 0.2"
)
#####################################################################################################################


########## parse command-line options ###############################################################################
newargs <- commandArgs(trailingOnly = TRUE)
if (length(newargs) == 0) {
  parser$parse_args("--help")
}

args <- parser$parse_args()

if (is.null(args$sub_name)) {
  stop("A subcommand is required: regression, survival, or complete", call. = FALSE)
}

require_nonempty <- function(value, option) {
  if (is.null(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(sprintf("Required option %s was not provided", option), call. = FALSE)
  }
}

validate_cli_file <- function(path, option) {
  require_nonempty(path, option)
  if (!file.exists(path)) {
    stop(sprintf("File supplied to %s does not exist: %s", option, path), call. = FALSE)
  }
}
#####################################################################################################################

########## REGRESSION ##########
if (args$sub_name == "regression") {
  validate_cli_file(args$fname, "-I/--filein")
  data <- read_reboot_table(args$fname, sep = "\t")

  rebootRegression(
    data = data,
    outprefix = args$out,
    bootstrap = args$boot_iter,
    groupsize = args$nel,
    percentagefilter = args$pf,
    type = args$type,
    variancefilter = args$var,
    force = args$force,
    seed = args$seed,
    ncores = args$ncores,
    followup = args$followup,
    plots = TRUE,
    table = TRUE,
    saveJSON = TRUE,
    saveRDS = TRUE,
    log = TRUE,
    report = FALSE
  )
}
################################

########## SURVIVAL ##########
if (args$sub_name == "survival") {
  validate_cli_file(args$fname, "-I/--filein")
  require_nonempty(args$sig, "-S/--signature")

  if (args$multi) {
    validate_cli_file(args$clin_file, "-C/--clinical")
  } else if (!is.null(args$clin_file) && nzchar(args$clin_file)) {
    stop("--clinical can only be used together with --multivariate", call. = FALSE)
  }

  data <- read_reboot_table(args$fname, sep = "\t")

  clindata <- if (args$multi) {
    read_reboot_table(args$clin_file, sep = "\t")
  } else {
    NULL
  }

  rebootSurvival(
    data = data,
    signature = args$sig,
    outprefix = args$out,
    multivariate = args$multi,
    clindata = clindata,
    roc = args$roc_curve,
    variancefilter = args$var,
    followup = args$followup,
    p.cutoff = args$p_cutoff,
    bootstrap = args$boot_iter,
    force = args$force,
    plots = TRUE,
    table = TRUE,
    saveJSON = TRUE,
    saveRDS = TRUE,
    log = TRUE,
    report = FALSE
  )
}
##############################

########## COMPLETE ##########
if (args$sub_name == "complete") {
  validate_cli_file(args$fname, "-I/--filein")

  if (args$multi) {
    validate_cli_file(args$clin_file, "-C/--clinical")
  } else if (!is.null(args$clin_file) && nzchar(args$clin_file)) {
    stop("--clinical can only be used together with --multivariate", call. = FALSE)
  }

  data <- read_reboot_table(args$fname, sep = "\t")

  clindata <- if (args$multi) {
    read_reboot_table(args$clin_file, sep = "\t")
  } else {
    NULL
  }

  rebootComplete(
    data = data,
    outprefix = args$out,
    bootstrap = args$boot_iter,
    groupsize = args$nel,
    percentagefilter = args$pf,
    variancefilter = args$var,
    type = args$type,
    multivariate = args$multi,
    clindata = clindata,
    roc = args$roc_curve,
    force = args$force,
    seed = args$seed,
    ncores = args$ncores,
    p.cutoff = args$p_cutoff,
    followup = args$followup,
    plots = TRUE,
    table = TRUE,
    saveJSON = TRUE,
    saveRDS = TRUE,
    log = TRUE,
    report = FALSE
  )
}
##############################
