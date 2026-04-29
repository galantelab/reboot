#!/usr/bin/env Rscript

#Pay carefully attention to help errors. "Argparse" apparently bugs depending on help text size

library(argparse)

#create the top-level parser
parser <- ArgumentParser()
subparsers <- parser$add_subparsers(dest = "sub_name", metavar = "<subcommand>", help = "choose only one option")

#version
parser$add_argument('-v', '--version', action = 'version', version = 'reboot 1.2.0')

########## create the parser for the "REGRESSION" command ##########
parser_reg <- subparsers$add_parser('regression', help = "generates a molecular signature through penalized LASSO Cox regression")

parser_reg$add_argument("-I", "--filein", type = "character", dest = "fname", metavar = '',
                        help = "Input file name. Tab separated values (tsv) file containing expression and survival parameters")
parser_reg$add_argument("-O", "--outprefix", type = "character", dest = "out", metavar = '', default = "reboot",
                        help = 'Output file prefix. Default: reboot')
parser_reg$add_argument("-B","--bootstrap", type = "integer", dest = "boot_iter",  metavar = '', default = "1",
                        help = "Number of iterations for bootstrap simulation (int). Default: 1")
parser_reg$add_argument("-G", "--groupsize", type = "integer", dest = "nel",  metavar = '', default = "10",
                        help = "Number of genes/transcripts for each bootstrap simulation (int). Default: 10")
parser_reg$add_argument("-P", "--pcentfilter", type = "numeric", dest = "pf",  metavar = '', default = "0.3",
                        help = "Percentage of correlated feature pairs allowed in each iteration (numeric). Default: 0.3")
parser_reg$add_argument("-T", "--type", type = "character", dest = "ty",  metavar = '', default = "gene",
                        help = 'Type of transcriptome data: gene or transcript (character). Default: gene')
parser_reg$add_argument("-V", "--varfilter", type = "numeric", dest = "var",  metavar = '', default = "0.01",
                        help = 'Minimum normalized variance (0-1) required for each feature among samples (numeric). Default: 0.01')
parser_reg$add_argument("-F", "--force", type = "logical", dest = "force", action = "store_true", default = FALSE,
                        help = "Choose -F to bypass OS and OS.time filters")
#####################################################################################################################

########## create the parser for the "SURVIVAL" command ##########
parser_sur <- subparsers$add_parser('survival', help = "applies a molecular signature in survival analysis for prognosis")

parser_sur$add_argument("-I", "--filein", type = "character", dest = "fname",  metavar = '',
                        help = "Input file name. Tab separated values (tsv) file containing expression and survival parameters")
parser_sur$add_argument("-O", "--outprefix", type = "character", dest = "out",  metavar = '', default = "reboot",
                        help = 'Output file prefix. Default: reboot')
parser_sur$add_argument("-S", "--signature", type = "character", dest = "sig",  metavar = '',
                        help = 'Tab separated values (tsv) file containing a set of features and corresponding regression coefficients')
parser_sur$add_argument("-M", "--multivariate", type = "logical", dest = "type", action = "store_true", default = FALSE,
                        help = 'If clinical variables should be included, choose -M. This option is tied with -C option')
parser_sur$add_argument("-C", "--clinical", type = "character", dest = "clin_file",  metavar = '', default = "",
                        help = 'Tab separated values (tsv) file containing binary categorical variables only. Required if -M option is chosen')
parser_sur$add_argument("-R", "--roc", type = "logical", dest = "roc_curve", action = "store_true", default = FALSE,
                        help = 'If molecular score should be categorized according to a ROC curve instead of the median, choose -R')
parser_sur$add_argument("-V", "--varfilter", type = "numeric", dest = "var",  metavar = '', default = "0.01",
                        help = 'Minimum normalized variance (0-1) required for follow up time (numeric). Default: 0.01')
parser_sur$add_argument("-F", "--force", type = "logical", dest = "force", action = "store_true", default = FALSE,
                        help = "Choose -F to bypass OS and OS.time filters")
#####################################################################################################################

########## create the parser for the "COMPLETE" command ##########
parser_all <- subparsers$add_parser('complete', help = "generates and applies a molecular signature (integrated analysis)")

parser_sur$add_argument("-I", "--filein", type = "character", dest = "fname",  metavar = '',
                        help = "Input file name. Tab separated values (tsv) file containing expression and survival parameters")
parser_sur$add_argument("-O", "--outprefix", type = "character", dest = "out",  metavar = '', default = "reboot",
                        help = 'Output file prefix. Default: reboot')
parser_sur$add_argument("-M", "--multivariate", type = "logical", dest = "type", action = "store_true", default = FALSE,
                        help = 'If clinical variables should be included, choose -M. This option is tied with -C option')
parser_sur$add_argument("-C", "--clinical", type = "character", dest = "clin_file",  metavar = '', default = "",
                        help = 'Tab separated values (tsv) file containing binary categorical variables only. Required if -M option is chosen')
parser_sur$add_argument("-R", "--roc", type = "logical", dest = "roc_curve", action = "store_true", default = FALSE,
                        help = 'If molecular score should be categorized according to a ROC curve instead of the median, choose -R')
parser_reg$add_argument("-B","--bootstrap", type = "integer", dest = "boot_iter",  metavar = '', default = "1",
                        help = "Number of iterations for bootstrap simulation (int). Default: 1")
parser_reg$add_argument("-G", "--groupsize", type = "integer", dest = "nel",  metavar = '', default = "10",
                        help = "Number of genes/transcripts for each bootstrap simulation (int). Default: 10")
parser_reg$add_argument("-P", "--pcentfilter", type = "numeric", dest = "pf",  metavar = '', default = "0.3",
                        help = "Percentage of correlated feature pairs allowed in each iteration (numeric). Default: 0.3")
parser_reg$add_argument("-T", "--type", type = "character", dest = "ty",  metavar = '', default = "gene",
                        help = 'Type of transcriptome data: gene or transcript (character). Default: gene')
parser_reg$add_argument("-V", "--varfilter", type = "numeric", dest = "var",  metavar = '', default = "0.01",
                        help = 'Minimum normalized variance (0-1) required for each feature and follow up time (numeric). Default: 0.01')
parser_reg$add_argument("-F", "--force", type = "logical", dest = "force", action = "store_true", default = FALSE,
                        help = "Choose -F to bypass OS and OS.time filters")
#####################################################################################################################

#Triggers help message if requested by user
args <- parser$parse_args()
newargs <- commandArgs(trailingOnly = TRUE)
if (length(newargs) == 0) {parser$parse_args('-h')}

########## REGRESSION ##########
if (args$sub_name == "regression") {
	system(paste(
		"regression.R",
		"-I", args$fname,
		"-O", args$out,
		"-B", args$boot_iter,
		"-G", args$nel,
		"-P", args$pf,
		"-T", args$ty,
		"-V", args$var,
		"-F", args$force,
		collapse = " "))
}
################################

########## SURVIVAL ##########
if (args$sub_name == "survival") {
	logfile <- paste(args$out, ".log", sep = "")
	if (file.exists(logfile)) {unlink(logfile)}

	#Check if clinical data is provided in case type == TRUE
	if (args$type) {
		if (args$clin_file == "") {
			cat("Insert file with clinical variables\n")
			q(status = 0)
		} else{
			system(paste(
			"survival.R",
			"-I", args$fname,
			"-O", args$out,
			"-S", args$sig,
			"-M", args$type,
			"-C", args$clin_file,
			"-R", args$roc_curve,
			"-V", args$var,
			"-F", args$force,
			collapse = " "))
		}
	} else {
		system(paste(
		"survival.R",		
		"-I", args$fname,
		"-O", args$out,
		"-S", args$sig,
		"-M", args$type,
		"-R", args$roc_curve,
		"-V", args$var,
		"-F", args$force,
		collapse = " "))
	}
}
##############################

########## COMPLETE ##########
if (args$sub_name == "complete") {

	#Check if clinical data is provided in case type == TRUE
	if (args$type) {
		if (args$clin_file == "") {
			cat("Insert file with clinical variables\n")
			q(status = 0)
		}
	}

	#Run "regression" - module I
  system(paste(
	"regression.R",
	"-I", args$fname,
	"-O", args$out,
	"-B", args$boot_iter,
	"-G", args$nel,
	"-P", args$pf,
	"-T", args$ty,
	"-V", args$var,
	"-F", args$force,
	collapse = " "),
	intern = TRUE)

  #Save signature automatically
  int_sig <- paste(args$out, "_signature.txt", sep = "")

  #Run "survival" - module II
  ##Check if clinical data is provided in case type == TRUE
	if (args$type) {
		if (args$clin_file == "") {
			cat("Insert file with clinical variables\n")
			q(status = 0)
		} else{
			system(paste(
			"survival.R",
			"-I", args$fname,
			"-O", args$out,
			"-S", int_sig,
			"-M", args$type,
			"-C", args$clin_file,
			"-R", args$roc_curve,
			"-V", args$var,
			"-F", args$force,
			collapse = " "))
		}
	} else {
		system(paste(
		"survival.R",
		"-I", args$fname,
		"-O", args$out,
		"-S", int_sig,
		"-M", args$type,
		"-R", args$roc_curve,
		"-V", args$var,
		"-F", args$force,
		collapse = " "))
	}
}
##############################
