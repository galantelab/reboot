#!/usr/bin/env Rscript

#suppressMessages(library(optparse))
#options(lifecycle_disable_verbose_retirement = TRUE)
#options(warn=-1)

#Time counter
start_time <- Sys.time()

######Argument parser######

option_list <- list(

optparse::make_option(c("-I", "--filein"), action="store",
type='character', dest = "fname", help="Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival paramenters"),

optparse::make_option(c("-O", "--outprefix"), action="store",
type='character', dest = "out", default = "reboot", help="Output file prefix. Default: reboot"),

optparse::make_option(c("-B","--bootstrap"), action = "store",
type = "integer", dest = "booty", default = "5", help = "Number of iterations for bootstrap simulation (int). Default: 5"),

optparse::make_option(c("-G", "--groupsize"), action="store",
type='integer', dest = "nel", default = "10", help="Number of genes/transcripts to be selected in each bootstrap simulation (int). Default: 10"),

optparse::make_option(c("-P", "--percentagefilter"), action="store",
type='numeric', dest = "pf", default = "0.3", help="Percentage of correlated gene/transcript pairs allowed in each iteration. Default: 0.3"),

optparse::make_option(c("-V", "--variancefilter"), action="store",
type='numeric', dest = "var", default = "0.01", help="Minimum normalized variance (0-1) required for each gene/transcript among samples (double). Default: 0.01"),

optparse::make_option(c("-T", "--type"), action="store",
type='character', dest = "ty", default = "gene", help="Declare which type of transcriptome data to be analyzed: gene or transcript. Default: gene"),

optparse::make_option(c("-F", "--force"), action="store",
type='logical', dest = "fierce", default = FALSE, help="To force overcome follow up variance filter and/or proportion filter for survival status (<20%), choose -F"))

opo <- optparse::OptionParser(option_list=option_list, add_help_option = T)
in_object <- optparse::parse_args(opo)
logname <- in_object$out
outname <- paste(in_object$out, "_signature.txt", sep="")
outplot <- in_object$out
fierce <-in_object$fierce
ty <- in_object$ty
 
####Importing libraries####

#suppressMessages(library("mice"))
#suppressMessages(library("penalized"))
#suppressMessages(library("tidyverse"))
#suppressMessages(library("hash"))
#suppressMessages(library("R.utils"))

######read input file######

full_data = read.table(in_object$fname, header=T, row.names=1, check.names=F)
colnames(full_data) <- gsub("-","__",colnames(full_data))
if (ty=="gene") {
	bar = 0.0035
} else if(ty=="transcript") {
	bar = 0.011
} else {
	cat("Type has to be either gene or transcript\n")
	q(status=0)
} 
####Setting plot theme#####

mytheme <- ggplot2::theme_bw() + 
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(face = "plain", colour = "black"),
        axis.text = ggplot2::element_text(face = "bold", colour = "black"),
        legend.text = ggplot2::element_text(colour = "black", face = "plain"),
        legend.title = ggplot2::element_text(colour = "black", face = "bold"),
        axis.ticks = ggplot2::element_line(colour = "black"), axis.line = ggplot2::element_line(colour = "black"))

#######Log file##########

sink(file = paste(logname, ".log", sep=''))

# To be added in the beginning of the script since this is intended to be the first lines of the ".log" file
cat("\n\n============================================================")
cat(" Make Signature ")
cat("============================================================\n\n")

scriptname <- "Regression analysis"
cat("Chosen parameters: ")
cat(paste(commandArgs(trailingOnly = T), collapse = " "))
cat("\n\n")

#####Checking data#########

nlines <- function(full_data, gs){
	if (nrow(full_data)<=10){
		sink(file = paste(logname, ".error", sep=''))
		cat("Error. There are less than 10 instances. Please, increase the number of lines for a proper analysis. \n")
		sink()
		q(status=0)
	}
	if (ncol(full_data) > 10) {
		if ((nrow(full_data) < 30) & (gs > 10)){
			cat("The proportion of instances per attributes might be low. You may consider to lower group size for a better analysis. \n\n")
		}

	}
}

#####Schoenfeld test######

ph_assumptions <- function(full_data){
	cat("Performing schoenfeld test\n\n")
	filt <- vector()
	attributes <- colnames(full_data[3:dim(full_data)[2]])
	for (i in attributes){
		phmodel <- survival::coxph(formula = formula(paste('survival::Surv(OS.time, OS)~', i)) , data = full_data)
		if(!is.na(phmodel$coef)){
			tryCatch({
				phmodel <- survival::coxph(formula = formula(paste('survival::Surv(OS.time, OS)~', i)) , data = full_data)
				schoen <- survival::cox.zph(phmodel)
				pval <- schoen$table[1,3]
				if (pval > 0.05){
					filt <- c(filt, i)
				}
			},warning=function(w){}, error=function(e){})
		}
	}
	losers <- setdiff(attributes, filt)
	cat(length(losers)," columns not allowed by schoenfeld test: ",losers, "\n\n")
	return(full_data[,c(colnames(full_data)[1:2],filt)])
	
}


######Error check 1#######

numberfilter1 <- function(dataf, g, outname, outplot) {
	if ((ncol(dataf)-2) < g){
		cat("The number of columns per group exceeds the number of columns", "\n", "\n")
		cat("Performing single multivariate regression","\n","\n")
		coemale <- regression(dataf)
		feature <- names(coemale)
		coefficient <- unname(coemale)
		if (any(!(coefficient==0))){
			coemale <- cbind(feature,coefficient)
			coemale <- gsub("__","-",coemale)
			#coemale <- gsub("@#!", "-", coemale)   #back to initial names
			write.table(coemale, outname, sep="\t", row.names=F, quote=F)
			#histogram(outplot,coemale)
			cat("Done", "\n")
		}
		else {
			cat("No signature found, all coefficients are equal 0", "\n")
			sink()
		}
		q(status=0)
	}
}

######Error check 2#######

numberfilter2 <- function(dataf, g, outname, outplot) {

	if (ncol(dataf) <= 2) {
		sink(file = paste(out, ".err", sep=''), append=T)
		cat("No column was left after variance filter", "\n","\n")
		sink()
		q(status=0)
	}

	if ((ncol(dataf) - 2) < g) {
		    cat("The number of columns is lower than group size due to variance filter","\n","\n")
		    cat("Performing single multivariate regression","\n","\n")
		    coemale <- regression(dataf)
		    feature <- names(coemale)
		    coefficient <- unname(coemale)
		    if (any(!(coefficient==0))){
                        coemale <- cbind(feature,coefficient)
			coemale <- gsub("__","-",coemale)
			#coemale <- gsub("@#!", "-", coemale)   #back to initial names
                        write.table(coemale, outname ,sep="\t", row.names=F, quote=F)
			histogram(outplot,coemale)
                        cat("Done", "\n")
                    }
                    else {
                        cat("No signature found, all coefficients are equal 0", "\n")
		    }
		    q(status=0)
       }
}
			
######Correlation filter######

corfun <- function(cmatrix, pf){
	indexes = c()
	pval = c()
	names= c()
	ngenes = 3:(ncol(cmatrix) - 1)
	for (t in ngenes){
		for (u in ((t + 1): (ncol(cmatrix)))){
			aux <- suppressWarnings(cor.test(x=cmatrix[,t], y=cmatrix[,u], method = 'spearman'))
			indexes <- c(indexes, round(aux$estimate,3))
		        pval <- c(pval, round(aux$p.value,3))
			names <- c(names, paste(colnames(cmatrix)[t],colnames(cmatrix)[u],sep="_"))
		}
	}
	if (((sum((indexes > 0.80) & (pval < 0.05))) / length(pval)) >= pf){ 
		switch=1
		cat("This iteration was avoided due to correlation among columns. Sperman correlation values and p-values are respectively:", "\n")
		cat(paste(names,indexes,pval,sep=":"),"\n","\n")

	}
	else{
		switch=0
	}
	return(switch)
}


######Bootfunction######

bootstrapfun <- function(full_data, booty, nel , outname, outplot, pf, bar){
	
	##setting up hash## 
	
	cat("Starting bootstrap ", booty, "iterations","\n\n")
	k <- colnames(full_data[3:length(colnames(full_data))])
	v <- vector("list", length(k))
	yield <- hash::hash(k,v)
	
	##looping##		
	i=1
	while (i<=booty){
		cat("processing iteration: ",i, "\n","\n")
		if(ncol(full_data) == 3) {
			cmatrix = full_data
		} else {
			cmatrix <- subsample(full_data, nel, i)
			#checking correlation#
			if (nel > 1){   #1 element avoided
				if (corfun(cmatrix, pf) == 1){
					next
				}
			}
		}
		#running regression#

		coemale <- regcall(cmatrix, nel, full_data) 
		##saving coeficients##	
	
		if (!is.null(coemale)){
			i = i + 1
			for (j in 1:length(coemale)){
				name = names(coemale[j])
				val = coemale[j][[1]]
				eval(parse(text=paste("yield$", name, "<- c(yield$", name, "," , val, ")", sep="")))
			}
		}
	}
  
	#Processing result
	aux=c(NULL,NULL, NULL)
	
	##calculating mean##

	for (feature in hash::keys(yield)){
	coefficient <- suppressWarnings(eval(parse(text=paste("mean(yield$", feature, ")", sep=""))))
        sd <- suppressWarnings(eval(parse(text=paste("sd(yield$", feature, ")", sep=""))))
	aux <- rbind(aux, cbind(feature, coefficient, sd))
	}	
	tt <- as.data.frame(aux)
	if (any(!complete.cases(tt$coefficient))){
		cat("NA coefficient found, increase coverage for a proper analysis", "\n")
	}

	tt <- dplyr::filter(tt, abs(as.numeric(as.character(coefficient))) >= bar)

	
	if (any(!(tt$coefficient == 0)) & dim(tt)[1]!=0){
		tt$feature <- gsub("__","-",tt$feature)
		write.table(tt, outname, sep="\t", row.names=F, quote=F)
		histogram(outplot,tt)
		lolli(outplot,tt)
	}
	else {
		cat("No signature found, all coefficients are not significant", "\n")
}
	cat("Done\n")

}

	
######Variance filter######

varfun <- function(cmatrix, var, file, fierce, out) {
	if (any(apply(cmatrix, 2, function(x) any(is.na(x))))){
		cat("Cheking NAs\n")
		impu <- mice::mice(cmatrix, print=F)
		cmatrix <-  mice::complete(impu)	
	}
	colnames(cmatrix)[1:2] <- c("OS","OS.time")	

	if (ncol(cmatrix) == 3) {
		divisor = max(cmatrix[3])
		dividendo = cmatrix[3]
		normalized = divisor/dividendo
		variance = var(normalized)
		if (variance < var) {
			cat("All columns rejected by variance filter","\n","\n")
	                q(status=0)
		}
		cat("No columns rejected by variance filter","\n","\n")
	} else {
	maxes <- matrix(apply(cmatrix[,3:ncol(cmatrix)],2,max), nrow=1)
	if (0 %in% maxes){
		cat("Columns with only 0s found in ", file, ". Remove such columns and try again.", "\n")
		q(status=0)
	}
	cat("Calculating normalized variances", "\n\n")
	dividendo <- dplyr::bind_rows(replicate(nrow(cmatrix), as.data.frame(maxes), simplify=F))
	divisor <- cmatrix[,3:ncol(cmatrix)]
	normalized <- divisor/dividendo
	variances <- apply(normalized,2,var)
	filtered <- c()
	losers <- c()
	for (i in (1:length(variances))) {
		if (variances[i] > var){
			filtered <- c(filtered,i)
		}
	}
	if (class(filtered)=="NULL"){
		cat("All columns rejected by variance filter","\n","\n")
		q(status=0)
	} else if(length(filtered) == length(variances)) {
		cat("No columns rejected by variance filter","\n","\n")	
	} else {
		losers <- names(variances[-filtered])
		filtered <- filtered+2
		cmatrix <- cmatrix[, c(1,2,filtered)]
		cat (length(losers)," columns with variance lower than ", var, " was removed from analysis: ",losers, "\n","\n")

	}
	#Dealing with SO and SO time
	if (!fierce){
		OSstatus <- cmatrix[,1]
		percentage <- sum(OSstatus)/length(OSstatus)
		if (percentage < 0.2 | percentage > 0.8){
			sink(file = paste(out, ".err", sep=''), append=T)
			cat("Survival status proportion:", percentage, " is probably not enough to the analysis. \nDeath or recidive are alternative options for the analysis", "\n\n")	
			cat("If you want to continue anyway, choose the flag F. \n\n")
			sink()
			q(status=0)
		}
		
		followup <- cmatrix[,2]
		uplimit <- max(followup)	
		normalized <- followup/uplimit
		fvar <- var(normalized)
		if (fvar < var){
			sink(file = paste(out, ".err", sep=''), append=T)
			cat("Follow up variance: ", var, " has not passed the variance test. \n\n")
			cat("If you want to continue anyway, choose the flag F, or lower variance filter. \n\n")
			sink()
			q(status=0)
		}
		
	}
	}
	return(cmatrix)  

}	


######Subsampling procedure#########

subsample <- function(full_data, nel, seed=i){
	set.seed(seed)
	shuffle <- sample(colnames(full_data[,3:ncol(full_data)]), size=nel, replace=F)
	cat("Picked columns: ",shuffle,"\n","\n")
	cmatrix <- cbind(full_data[,1:2],subset(full_data,select=shuffle))
		
}


######Regression time keeper#######

regcall <- function(cmatrix, nel, full_data){
	tryCatch(R.utils::withTimeout(coemale <- regression(cmatrix), timeout=60, onTimeout= "warning"),
		warning=function(warning_condition){
			cat("Regression time exceeded, you may consider changing variance and/or correlation filters. Trying again \n");
			cmatrix <- subsample(full_data,nel);
			regcall(cmatrix, nel, full_data)	
		},
		error=function(e){
			#coemale <- NULL
			#return(coemale)
		    }
	)
}	


######Regression procedure############

regression <- function(cmatrix){

	#fit1 <- penalized::profL1(survival::Surv(OS.time,OS)~., data=cmatrix, fold=10, maxlambda1=100, plot=F, trace=F)
	#options(show.error.messages = F)
        try(fit1 <- penalized::profL1(survival::Surv(OS.time,OS)~., data=cmatrix, fold=10, plot=F, trace=F))
	fit1 <- penalized::profL1(survival::Surv(OS.time,OS)~., data=cmatrix, fold=10, plot=F, trace=F)
	#fit2 <- penalized::profL2(survival::Surv(OS.time,OS)~., data=cmatrix, fold=fit1$fold, minl = 0.1, maxlambda2 = 10)
	#opt1 <- penalized::optL1(survival::Surv(OS.time,OS)~., data=cmatrix, fold=fit1$fold, maxlambda1=10, trace=F )
	opt1 <- penalized::optL1(survival::Surv(OS.time,OS)~., data=cmatrix, fold=fit1$fold, trace=F )
	#opt2 <- penalized::optL2(survival::Surv(OS.time,OS)~., data=cmatrix, fold=fit2$fold)
	fit <- penalized::penalized(survival::Surv(OS.time,OS)~., data=cmatrix, lambda1=opt1$lambda, trace=F)
	coemale <- penalized::coefficients(fit, "all")
	return(coemale)
	
}


#######lollipop plot#################

lolli <- function(out,tt){
	cat("building ranking plot\n")
	#fname <- paste(out,"_lollipop.png",sep="")
	fname <- paste(out,"_lollipop.pdf",sep="")
	tt<-tt[complete.cases(tt), ]
	tt$coefficient <- as.numeric(as.character(tt$coefficient))
	tt <- dplyr::filter(tt,coefficient!=0)
	tt$feature = factor(tt$feature, levels=tt[order(tt$coefficient),"feature"])
	filter <- as.character(dplyr::top_n(tt, 10, abs(coefficient))$feature)
	tt <- dplyr::filter(tt, feature %in% filter)
	tt <- dplyr::mutate(tt, name = forcats::fct_reorder(feature, coefficient))

	pdf(fname)					
		ll <- ggplot2::ggplot(tt, ggplot2::aes(x=reorder(feature, coefficient), y=coefficient)) +
		ggplot2::geom_segment(ggplot2::aes(x=feature, xend=feature, y=0, yend=coefficient), size=2, color="grey") +
		ggplot2::geom_point(size=4, colour="#1c9099") +
		#ggplot2::theme_light(base_family = "Helvetica") +
		mytheme +
		ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
			axis.ticks.x = ggplot2::element_blank(),
			axis.title.y = ggplot2::element_blank(),
			axis.title.x = ggplot2::element_blank(),
			panel.border = ggplot2::element_blank(),
			#axis.text.y = ggplot2::element_text(angle=90),
			#axis.text.x = ggplot2::element_text(face="bold",angle=45, size=12, hjust=0),
			axis.text.x = ggplot2::element_text(angle=60, hjust=1))
			#axis.text.y= ggplot2::element_text(face="bold")) 
	#scale_y_continuous(breaks = round(seq(min(tt$coefficient), max(tt$coefficient), by = 0.005),3)) 
	print (ll, newpage=F)
	dev.off()
	
	
}

#######Histogram plot#################

histogram <- function(out,tt){
	cat("building histogram plot\n")
	fname<-paste(out,"_hist.pdf",sep="")
	tt<-tt[complete.cases(tt), ]
	tt$coefficient <- as.numeric(as.character(tt$coefficient))
	tt <- dplyr::filter(tt, coefficient != 0)
	pdf(fname)
	pl<-ggplot2::ggplot(tt, ggplot2::aes(x = coefficient))+
		ggplot2::geom_histogram(binwidth = 0.005, alpha=1, position="identity") +
		ggplot2::scale_y_continuous(expand = c(0,0)) +
		ggplot2::xlab("coefficients")+
		ggplot2::ylab("") +
		mytheme
		print (pl)
	dev.off()
}

####Main####

#check number of lines#

#nlines(full_data, in_object$nel)

#Perform variance filter#

full_data <- varfun(full_data, in_object$var, in_object$fname, fierce, logname)

#Perform schoenfeld tests#

full_data <- ph_assumptions(full_data)

#check file error#

numberfilter1(full_data, in_object$nel, outname, outplot)

#check file error 2#

numberfilter2(full_data, in_object$nel, outname, outplot)

#Perform Regression#

bootstrapfun(full_data, in_object$booty, in_object$nel, outname, outplot, in_object$pf, bar)


####Time feedback####

end_time <- Sys.time()
elapsed_time <- difftime(time1 = end_time, time2 = start_time, units = "secs")

if (elapsed_time >= 3600) {
  cat(paste("\nTime to run ", scriptname, ": ", round(x = (elapsed_time[[1]] / 3600), digits = 2), " hours\n", sep = ""))
} else {
  if (elapsed_time >= 60) {
    cat(paste("\nTime to run ", scriptname, ": ", round(x = (elapsed_time[[1]] / 60), digits = 2), " minutes\n", sep = ""))
  } else {
    cat(paste("\nTime to run ", scriptname, ": ", round(x = elapsed_time[[1]], digits = 2), " seconds\n", sep = ""))
  }
}

sink()
