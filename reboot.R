#!/usr/local/bin/Rscript

library(argparse)

#create the top-level parser
parser = ArgumentParser()
#parser$add_argument()
subparsers = parser$add_subparsers(dest="sub_name", metavar = "<subcommand>", help = "choose only one option")

#version
parser$add_argument('-v', '--version', action='version', version='reboot 1.0.0')


#create the parser for the "regression" command

parser_reg = subparsers$add_parser('regression', help= "generates signature through multivariate cox regression")

parser_reg$add_argument("-I", "--filein", type="character", dest = "fname", metavar ='',
			help="Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival data")

parser_reg$add_argument("-O", "--outprefix", type="character", dest = "out", metavar = '',
			default = "reboot", help='Output file prefix (str). Default: reboot')

parser_reg$add_argument("-B","--bootstrap", type = "integer", dest = "booty",  metavar = '',
                        default = "1", 
                        help = "Number of iterations for bootstrap simulation (int). Default: 1")

parser_reg$add_argument("-G", "--groupsize", type="integer", dest = "nel",  metavar = '',
                        default="10",
                        help="Number of genes/transcripts to be selected in each bootstrap simulation (int). Default: 10")

parser_reg$add_argument("-P", "--pcentfilter", type="character", dest = "pf",  metavar = '',
                        default = "0.3",
                        help="Percentage of correlated gene/transcript pairs allowed in each iteration (double). Default: 0.3")

parser_reg$add_argument("-V", "--varfilter", type="character", dest = "var",  metavar = '',
                        default = "0.01",
                        help="Minimum normalized variance (0-1) required for each gene/transcript among samples and follow up time (double). Default: 0.01")

parser_reg$add_argument("-F", "--force", 
                        dest = "force", action="store_true",
                        default = FALSE,
                        help="To force overcome follow up variance filter and/or proportion filter for survival status (<20%), choose -F")

#create the parser for the "survival" command

parser_sur = subparsers$add_parser('survival', help = "applies sgnature in survival analysis")


parser_sur$add_argument("-I", "--filein", type="character", dest = "fname",  metavar = '',
                       help="Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival paramenters")      

parser_sur$add_argument("-O", "--outprefix", type="character", dest = "out",  metavar = '',
                        default = "reboot",
                        help='Output file prefix. Default: reboot')

parser_sur$add_argument("-M", "--multivariate", 
                        dest = "type", action="store_true", default="FALSE",
                        help='If clinical variables should be included, choose -M. This option is tied with -C option')

parser_sur$add_argument("-C", "--clinical", type="character", dest = "clin_file",  metavar = '',
			default="",
                        help='Tab separated values (tsv) file containing binary categorical variables only. Required if -M option is chosen')

parser_sur$add_argument("-R", "--roc",
                        dest = "roc_curve", action="store_true", default="FALSE",
                        help='If genetic score should be categorized according to a ROC curve instead of the median, choose -R')

parser_sur$add_argument("-S", "--signature",  metavar = '',
                        type="character", dest = "sig",
                        help='Tab separated values (tsv) file containing a set of genes/transcripts and corresponding cox coefficients')

parser_all$add_argument("-V", "--varfilter", type="character", dest = "var",  metavar = '',
                        default = "0.01",
                        help="Minimum normalized variance (0-1) required for follow up time (double). Default: 0.01")

parser_sur$add_argument("-F", "--force", 
			dest = "force", action="store_true",
                        default = FALSE,
                        help="To force overcome follow up variance filter and/or proportion filter for survival status (<20%), choose -F")

#create the parser for the "complete" command

parser_all = subparsers$add_parser('complete', help = "generates and applies signature (integrated analysis)")

parser_all$add_argument("-I", "--filein", type="character", dest = "fname", metavar ='',
                        help='Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival data')

parser_all$add_argument("-O", "--outprefix", type="character", dest = "out", metavar = '',
                        default = "reboot", help='Output file prefix. Default: reboot')

parser_all$add_argument("-B","--bootstrap", type = "integer", dest = "booty",  metavar = '',
                        default = "1",
                        help = "Number of iterations for bootstrap simulation (int). Default: 1")

parser_all$add_argument("-G", "--groupsize", type="integer", dest = "nel",  metavar = '',
                        default="10",
                        help="Number of genes/transcripts to be selected in each bootstrap simulation (int). Default: 10")

parser_all$add_argument("-P", "--pcentfilter", type="character", dest = "pf",  metavar = '',
                        default = "0.3",
                        help="Percentage of correlated gene/transcript pairs allowed in each iteration. Default: 0.3")

parser_all$add_argument("-V", "--varfilter", type="character", dest = "var",  metavar = '',
                        default = "0.01",
                        help="Minimum normalized variance (0-1) required for each gene/transcript among samples and follow up time (double). Default: 0.01")

parser_all$add_argument("-M", "--multivariate",
                        dest = "type", action="store_true", default=FALSE,
                        help='If clinical variables should be included, choose -M. This option is tied with -C option')

parser_all$add_argument("-C", "--clinical", type="character", dest = "clin_file",  metavar = '',
                        default="",
                        help='Tab separated values (tsv) file containing binary categorical variables only. Required if -M option is chosen')

parser_all$add_argument("-R", "--roc",
                        dest = "roc_curve", action="store_true", default=FALSE,
                        help='If continuous variables should be categorized according to a ROC curve instead of median, choose -R')


parser_all$add_argument("-F", "--force", 
                        dest = "force", action="store_true",
                        default = FALSE,
                        help="To force overcome follow up variance filter and/or proportion filter for survival status (<20%), choose -F")

#parse some argument lists
args = parser$parse_args()

#return help in case of empty call
newargs = commandArgs(trailingOnly=TRUE)
if (length(newargs)==0){
	parser$parse_args('-h')
}


if (args$sub_name=="regression"){
	system(paste(
#		"Rscript regression.R",
		"regression.R",
		"-I", args$fname,
		"-O", args$out,
		"-B", args$booty,
		"-G", args$nel,
		"-P", args$pf,
		"-V", args$var,
		"-f", arg$force,
		collapse=" "))
}

if (args$sub_name=="survival"){
	logfile <- paste(args$out,".log",sep="")
	if(file.exists(logfile)){
		unlink(logfile)
	}	


	#Check if clinical data is provided in case type == TRUE
	if (args$type){
		if (args$clin_file == ""){
			cat("Insert file with clinical variables\n")
			q(status=0)
		}else{		
			system(paste(
#			"Rscript survival.R",
			"survival.R",
			"-I", args$fname,
			"-O", args$out,
			"-S", args$sig,
			"-M", args$type,
			"-C", args$clin_file,
			"-R", args$roc_curve,
			"-f", arg$force,
			collapse=" "))
		}
	}else{
		system(paste(
#		"Rscript survival.R",
		"survival.R",		
		"-I", args$fname,
		"-O", args$out,
		"-S", args$sig,
		"-M", args$type,
		"-R", args$roc_curve,
		"-f", arg$force,
		collapse=" "))
	
	}
}

if (args$sub_name=="complete"){

	#Check if clinical data is provided in case type == TRUE
	if(args$type){
		if(args$clin_file == ""){
			cat("Insert file with clinical variables\n")
			q(status=0)
		}
	}	
	
	system(paste(
#	"Rscript regression.R",
	"regression.R",
	"-I", args$fname,
	"-O", args$out,
	"-B", args$booty,
	"-G", args$nel,
	"-P", args$pf,
	"-V", args$var,
	"-f", arg$force,
	collapse=" "),
	intern=TRUE)
	
	assinatura <- paste(args$out,"_signature.txt",sep="")
	#Check if clinical data is provided in case type == TRUE
	if (args$type){
		if (args$clin_file == ""){
			cat("Insert file with clinical variables\n")
			q(status=0)
		}else{		
			system(paste(
#			"Rscript survival.R",
			"survival.R",
			"-I", args$fname,
			"-O", args$out,
			"-S", assinatura,
			"-M", args$type,
			"-C", args$clin_file,
			"-R", args$roc_curve,
			"-f", arg$force,
			collapse=" "))
		}
	}else{
		system(paste(
#		"Rscript survival.R",
		"survival.R",
		"-I", args$fname,
		"-O", args$out,
		"-S", assinatura,
		"-M", args$type,
		"-R", args$roc_curve,
		"-f", arg$force,
		collapse=" "))
	
	}
}



