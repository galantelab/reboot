#!/usr/local/bin/Rscript

suppressMessages(library(optparse))
#options(lifecycle_disable_verbose_retirement = TRUE)
#options(warn=-1)

#Time counter
start_time <- Sys.time()

######Argument parser######

option_list <- list(

make_option(c("-I", "--filein"), action="store",
type='character', dest = "fname", help="Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival paramenters"),

make_option(c("-O", "--outprefix"), action="store",
type='character', dest = "out", default = "reboot", help="Output file prefix. Default: reboot"),

make_option(c("-B","--bootstrap"), action = "store",
type = "integer", dest = "booty", default = "5", help = "Number of iterations for bootstrap simulation (int). Default: 5"),

make_option(c("-G", "--groupsize"), action="store",
type='integer', dest = "nel", default = "10", help="Number of genes/transcripts to be selected in each bootstrap simulation (int). Default: 10"),

make_option(c("-P", "--percentagefilter"), action="store",
type='numeric', dest = "pf", default = "0.3", help="Percentage of correlated gene/transcript pairs allowed in each iteration. Default: 0.3"),

make_option(c("-V", "--variancefilter"), action="store",
type='numeric', dest = "var", default = "0.01", help="Minimum normalized variance (0-1) required for each gene/transcript among samples (double). Default: 0.01"))


opo <- OptionParser(option_list=option_list, add_help_option = T)
in_object <- parse_args(opo)
logname <- in_object$out
outname <- paste(in_object$out, "_signature.txt", sep="")
outplot <- in_object$out

#####Importing libraries####

suppressMessages(library("penalized"))
suppressMessages(library("tidyverse"))
suppressMessages(library("hash"))
suppressMessages(library("R.utils"))

######read input file######

full_data = read.table(in_object$fname, header=T, row.names=1, check.names=F)
colnames(full_data) <- gsub("-","__",colnames(full_data))


####Setting plot theme#####

mytheme <- theme_bw() + 
  theme(panel.grid.major.x = element_blank(), text = element_text(face = "plain", colour = "black"),
        axis.text = element_text(face = "bold", colour = "black"),
        legend.text = element_text(colour = "black", face = "plain"),
        legend.title = element_text(colour = "black", face = "bold"),
        axis.ticks = element_line(colour = "black"), axis.line = element_line(colour = "black"))



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

#####Schoenfeld test######

ph_assumptions <- function(full_data){
	cat("Performing schoenfeld test\n\n")
	filt <- vector()
	attributes <- colnames(full_data[3:dim(full_data)[2]])
	for (i in attributes){
		phmodel <- coxph(formula = formula(paste('Surv(OS.time, OS)~', i)) , data = full_data)
		schoen <- cox.zph(phmodel)
		pval <- schoen$table[1,3]
		if (pval > 0.05){
			filt <- c(filt, i)
		}
	} 
	losers <- setdiff(attributes, filt)
	cat(length(losers)," columns not allowed by schoenfeld test: ",losers, "\n\n")
	return(full_data)
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
			#coemale <- gsub("@#!", "-", coemale)   #back to initial names
			write.table(coemale, outname, sep="\t", row.names=F, quote=F)
			#histogram(outplot,coemale)
			cat("Done", "\n")
		}
		else {
			cat("No signature found, all coefficients are equal 0", "\n")
		}
		q(status=0)
	}
}

######Error check 2#######

numberfilter2 <- function(dataf, g, outname, outplot) {

	if (ncol(dataf) <= 2) {
		cat("No column was left after variance filter", "\n","\n")
		q(status=0)
	}

	if ((ncol(dataf)-2) < g) {
		    cat("The number of columns is lower than group size due to variance filter","\n","\n")
		    cat("Performing single multivariate regression","\n","\n")
		    coemale <- regression(dataf)
		    feature <- names(coemale)
		    coefficient <- unname(coemale)
		    if (any(!(coefficient==0))){
                        coemale <- cbind(feature,coefficient)
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

corfun <- function(male_data, pf){
	indexes = c()
	pval = c()
	names= c()
	ngenes = 3:(ncol(male_data)-1)
	for (t in ngenes){
		for (u in ((t+1): (ncol(male_data)))){
			aux <- suppressWarnings(cor.test(x=male_data[,t], y=male_data[,u], method = 'spearman'))
			indexes <- c(indexes, round(aux$estimate,3))
		        pval <- c(pval, round(aux$p.value,3))
			names <- c(names, paste(colnames(male_data)[t],colnames(male_data)[u],sep="_"))
		}
	}
	if (((sum((indexes>0.80) & (pval<0.05)))/length(pf))>=pf){ 
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

bootstrapfun <- function(full_data, booty, nel , outname, outplot, pf){
	
	##setting up hash## 
	
	cat("Starting bootstrap ", booty, "iterations","\n\n")
	k <- colnames(full_data[3:length(colnames(full_data))])
	v <- vector("list", length(k))
	yield <- hash(k,v)
	
	##looping##		
	i=1
	while (i<=booty){
	cat("processing iteration: ",i, "\n","\n")
	male_data <- subsample(full_data, nel)
	
	#checking correlation#
	if (corfun(male_data, pf)==1){
		next
	}
	
	i=i+1

	#running regression#
	
	coemale <- regcall(male_data, nel, full_data) 
	
	##saving coeficients##	

	for (j in 1:length(coemale)){
	name = names(coemale[j])
	val = coemale[j][[1]]
	eval(parse(text=paste("yield$", name, "<- c(yield$", name, "," , val, ")", sep="")))
	} 
	}  
	#Processing result
	aux=c(NULL,NULL)
	
	##calculating mean##

	for (feature in keys(yield)){
	coefficient <- suppressWarnings(eval(parse(text=paste("mean(yield$", feature, ")", sep=""))))
	aux <- rbind(aux, cbind(feature, coefficient))
	}	
	tt <- as.data.frame(aux)

	tt <- tt %>%
		filter(coefficient!=0)

	if (any(!complete.cases(tt$coefficient))){
		cat("NA coefficient found, increase coverage for a proper analysis", "\n")
	}

	if (any(!(tt==0))){
		tt$feature <- gsub("__","-",tt$feature)
		write.table(tt, outname, sep="\t", row.names=F, quote=F)
		histogram(outplot,tt)
		lolli(outplot,tt)
	}
	else {
		cat("No signature found, all coefficients are equal 0", "\n")
	}
	cat("Done\n")

}

	
######Variance filter######

varfun <- function(male_data, var, file) {
	maxes <- matrix(apply(male_data[,3:ncol(male_data)],2,max), nrow=1)
	if (0 %in% maxes){
		cat("Columns with only 0s found in ", file, ". Remove such columns and try again.", "\n")
		q(status=0)
	}
	cat("Calculating normalized variances", "\n\n")
	dividendo <- bind_rows(replicate(nrow(male_data), as.data.frame(maxes), simplify=F))
	divisor <- male_data[,3:ncol(male_data)]
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
		male_data <- male_data[, c(1,2,filtered)]
		cat (length(losers)," columns with variance lower than ", var, " was removed from analysis: ",losers, "\n","\n")

	}
	return(male_data)  

}	


######Subsampling procedure#########

subsample <- function(full_data, nel){

	shuffle <- sample(colnames(full_data[,3:ncol(full_data)]), size=nel, replace=F)
	cat("Picked columns: ",shuffle,"\n","\n")
	male_data <- cbind(full_data[,1:2],subset(full_data,select=shuffle))
		
}


######Regression time keeper#######

regcall <- function(male_data, nel, full_data){
	tryCatch(withTimeout(coemale <- regression(male_data), timeout=60, onTimeout= "warning"),
		warning=function(warning_condition){
			cat("Regression time exceeded, you may consider changing variance and/or correlation filters. Trying again \n");
			male_data <- subsample(full_data,nel);
			regcall(male_data, nel, full_data)	
		}
	)
}	


######Regression procedure############

regression <- function(male_data){

	#fit1 <- profL1(Surv(OS.time,OS)~., data=male_data, fold=10, maxlambda1=100, plot=F, trace=F)
	fit1 <- profL1(Surv(OS.time,OS)~., data=male_data, fold=10, plot=F, trace=F)
	#fit2 <- profL2(Surv(OS.time,OS)~., data=male_data, fold=fit1$fold, minl = 0.1, maxlambda2 = 10)
	#opt1 <- optL1(Surv(OS.time,OS)~., data=male_data, fold=fit1$fold, maxlambda1=10, trace=F )
	opt1 <- optL1(Surv(OS.time,OS)~., data=male_data, fold=fit1$fold, trace=F )
	#opt2 <- optL2(Surv(OS.time,OS)~., data=male_data, fold=fit2$fold)
	fit <- penalized(Surv(OS.time,OS)~., data=male_data, lambda1=opt1$lambda, trace=F)
	coemale <- coefficients(fit, "all")
	return(coemale)
	
}


#######lollipop plot#################

lolli <- function(out,tt){
	cat("building ranking plot\n")
	#fname <- paste(out,"_lollipop.png",sep="")
	fname <- paste(out,"_lollipop.pdf",sep="")
	tt<-tt[complete.cases(tt), ]
	tt$coefficient <- as.numeric(as.character(tt$coefficient))
	tt <- filter(tt,coefficient!=0)
	tt$feature = factor(tt$feature, levels=tt[order(tt$coefficient),"feature"])
	filter <- as.character(top_n(tt, 10, abs(coefficient))$feature)
	tt <- filter(tt, feature %in% filter)
	tt <- mutate(tt, name = fct_reorder(feature, coefficient))

	pdf(fname)					
		ll <- ggplot(tt, aes(x=reorder(feature, coefficient), y=coefficient)) +
		geom_segment(aes(x=feature, xend=feature, y=0, yend=coefficient), size=2, color="grey") +
		geom_point(size=4, colour="#1c9099") +
		#theme_light(base_family = "Helvetica") +
		mytheme +
		theme(panel.grid.major.x = element_blank(),
			axis.ticks.x = element_blank(),
			axis.title.y = element_blank(),
			axis.title.x = element_blank(),
			panel.border = element_blank(),
			#axis.text.y = element_text(angle=90),
			#axis.text.x = element_text(face="bold",angle=45, size=12, hjust=0),
			axis.text.x = element_text(angle=60, hjust=1))
			#axis.text.y= element_text(face="bold")) 
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
	tt <- filter(tt, coefficient != 0)
	pdf(fname)
	pl<-ggplot(tt, aes(x = coefficient))+
		geom_histogram(binwidth = 0.005, alpha=1, position="identity") +
		scale_y_continuous(expand = c(0,0)) +
		xlab("coefficients")+
		ylab("") +
		mytheme
		print (pl)
	dev.off()
}

####Main####

#Perform variance filter#

full_data <- varfun(full_data, in_object$var, in_object$fname)

#Perform schoenfeld tests#

full_data <- ph_assumptions(full_data)

#check file error#

numberfilter1(full_data, in_object$nel, outname, outplot)

#check file error 2#

numberfilter2(full_data, in_object$nel, outname, outplot)

#Perform Regression#

bootstrapfun(full_data, in_object$booty, in_object$nel, outname, outplot, in_object$pf)


####Time feedback####

end_time <- Sys.time()
elapsed_time <- difftime(time1 = end_time, time2 = start_time, units = "secs")

if (elapsed_time >= 3600) {
  cat(paste("Time to run ", scriptname, ": ", round(x = (elapsed_time[[1]] / 3600), digits = 2), " hours\n", sep = ""))
} else {
  if (elapsed_time >= 60) {
    cat(paste("Time to run ", scriptname, ": ", round(x = (elapsed_time[[1]] / 60), digits = 2), " minutes\n", sep = ""))
  } else {
    cat(paste("Time to run ", scriptname, ": ", round(x = elapsed_time[[1]], digits = 2), " seconds\n", sep = ""))
  }
}

sink()
