#' Lasso regression package for gene/transcript analysis 
#'
#' This function allows you to obtain beta coefficients of significance for genes or transcripts.
#' @param filein Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival paramenters.
#' @param outprefix Output file prefix. Default: reboot.
#' @param bootstrap Number of iterations for bootstrap simulation (int).Default: 1
#' @param groupsize Number of genes/transcripts to be selected in each bootstrap simulation (int). Default: 3
#' @param pcentfilter Percentage of correlated gene/transcript pairs allowed in each iteration. Default: 0.3
#' @param varfilter Minimum normalized variance (0-1) required for each gene/transcript among samples (double). Default: 0.01
#' @keywords regression genes transcripts
#' @export
#' @examples

get_coefficients <- function(filein, 
			     outprefix="reboot", 
			     bootstrap=1,
			     groupsize=3, 
			     pcentfilter=0.3, 
			     varfilter=0.01) {

	#Time counter
	start_time <- Sys.time()

	######read input file######

	full_data = read.table(filein, header=T, row.names=1, check.names=F)
	colnames(full_data) <- gsub("-","__",colnames(full_data))


	####Setting plot theme#####

	mytheme <- ggplot2::theme_bw() + 
	ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(face = "plain", colour = "black"),
		axis.text = ggplot2::element_text(face = "bold", colour = "black"),
	        legend.text = ggplot2::element_text(colour = "black", face = "plain"),
		legend.title = ggplot2::element_text(colour = "black", face = "bold"),
	        axis.ticks = ggplot2::element_line(colour = "black"), axis.line = ggplot2::element_line(colour = "black"))

	########Log file##########

	sink(file = paste(outprefix, ".log", sep=''))

	# To be added in the beginning of the script since this is intended to be the first lines of the ".log" file
	cat("\n\n============================================================")
	cat(" Make Signature ")
	cat("============================================================\n\n")

	scriptname <- "Regression analysis"
	cat("Chosen parameters: ")
	cat(paste(commandArgs(trailingOnly = T), collapse = " "))
	cat("\n\n")


	#Main#

	#check file error#
	
	numberfilter1(full_data, groupsize, outprefix, mytheme)
	
	#Perform variance filter#
	
	full_data <- varfun(full_data, varfilter, filein)
	
	#check file error 2#
	
	numberfilter2(full_data, groupsize, outprefix, mytheme)
	
	#Perform Regression#
	
	bootstrapfun(full_data, bootstrap, groupsize, outprefix, pcentfilter, mytheme)
	

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

}


######Error check 1#######

numberfilter1 <- function(dataf, g, out, mytheme) {
	if ((ncol(dataf)-2) < g){
		cat("The number of columns per group exceeds the number of columns", "\n", "\n")
		cat("Performing single multivariate regression","\n","\n")
		coemale <- regression(dataf)
		feature <- names(coemale)
		coefficient <- unname(coemale)
		if (any(!(coefficient==0))){
			coemale <- cbind(feature,coefficient)
			#coemale <- gsub("@#!", "-", coemale)   #back to initial names
			write.table(coemale, paste(out, "signature.txt", sep="_"), sep="\t", row.names=F, quote=F)
			histogram(out,coemale, mytheme)
			cat("Done", "\n")
		}
		else {
			cat("No signature found, all coefficients are equal 0", "\n")
		}
		q(status=0)
	}
}

######Error check 2#######

numberfilter2 <- function(dataf, g, out, mytheme) {

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
                        write.table(coemale, paste(out, "signature.txt", sep="_") ,sep="\t", row.names=F, quote=F)
			histogram(out,coemale, mytheme)
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

bootstrapfun <- function(full_data, booty, nel , out, pf, mytheme){
	
	##setting up hash## 
	
	cat("Starting bootstrap ", booty, "iterations","\n")
	k <- colnames(full_data[3:length(colnames(full_data))])
	v <- vector("list", length(k))
	yield <- hash::hash(k,v)
	
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
		
		#saving coeficients#	
		for (j in 1:length(coemale)){
			name = names(coemale[j])
			val = coemale[j][[1]]
			eval(parse(text=paste("yield$", name, "<- c(yield$", name, "," , val, ")", sep="")))
		} 
		print("teste4")
	}  
	
	#Processing result
	aux=c(NULL,NULL)
	
	##calculating mean##

	for (feature in hash::keys(yield)){
		coefficient <- suppressWarnings(eval(parse(text=paste("mean(yield$", feature, ")", sep=""))))
		aux <- rbind(aux, cbind(feature, coefficient))
	}	
	tt <- as.data.frame(aux)
	
	if (any(!(tt==0))){
		tt$feature <- gsub("__","-",tt$feature)
		write.table(tt, paste(out, "signature.txt", sep="_"), sep="\t", row.names=F, quote=F)
		histogram(out,tt, mytheme)
		lolli(out,tt, mytheme)
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
	dividendo <- dplyr::bind_rows(replicate(nrow(male_data), as.data.frame(maxes), simplify=F))
	divisor <- male_data[,3:ncol(male_data)]
	normalized <- divisor/dividendo
	variances <- apply(normalized,2,var)
	trash <- c()
	for (i in (1:length(variances))) {
		variances[i]
		if (variances[i] < var) {
			trash <- c(trash,names(variances[i]))
		}
	}

  if (class(trash)!="NULL"){
  male_data <- dplyr::select(male_data, -trash)
    cat (length(trash)," columns with variance lower than: ", var, " was removed from analysis: ",trash, "\n","\n")
  }
  else {
    cat("No columns rejected by variance filter","\n","\n")	
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
	print("teste1")
	tryCatch(R.utils::withTimeout(coemale <- regression(male_data), timeout=60, onTimeout= "warning"),
		warning=function(warning_condition){
			cat("Regression time exceeded, you may consider changing variance and/or correlation filters. Trying again \n");
			male_data <- subsample(full_data,nel);
			regcall(male_data, nel, full_data)	
		}
	)
}	


######Regression procedure############

regression <- function(male_data){
	print ("teste2")
	fit1 <- penalized::profL1(survival::Surv(OS.time,OS)~., data=male_data, fold=10, maxlambda1=100, plot=F, trace=F)
	#fit2 <- penalized::profL2(survival::Surv(OS.time,OS)~., data=male_data, fold=fit1$fold, minl = 0.1, maxlambda2 = 10)
	opt1 <- penalized::optL1(survival::Surv(OS.time,OS)~., data=male_data, fold=fit1$fold, maxlambda1=10, trace=F )
	#opt2 <- penalized::optL2(survival::Surv(OS.time,OS)~., data=male_data, fold=fit2$fold)
	fit <- penalized::penalized(survival::Surv(OS.time,OS)~., data=male_data, lambda1=opt1$lambda, trace=F)
	print("teste3")
	coemale <- penalized::coefficients(fit, "all")
	print("teste4")
	return(coemale)
	
}


#######lollipop plot#################

lolli <- function(out,tt, mytheme){
	cat("building ranking plot\n")
	#fname <- paste(out,"_lollipop.png",sep="")
	fname <- paste(out,"_lollipop.pdf",sep="")
	tt <- tt[complete.cases(tt), ]
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
		#theme_light(base_family = "Helvetica") +
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

histogram <- function(out,tt, mytheme){
	cat("building histogram plot\n")
	fname <- paste(out,"_hist.pdf",sep="")
	tt <- tt[complete.cases(tt), ]
	tt <- dplyr::filter(tt,coefficient!=0)
	tt$coefficient <- as.numeric(as.character(tt$coefficient))
	pdf(fname)
	pl <- ggplot2::ggplot(tt, ggplot2::aes(x=coefficient))+
		ggplot2::geom_histogram(binwidth = 0.005, alpha=1, position="identity") +
		ggplot2::scale_y_continuous(expand=c(0,0)) +
		ggplot2::xlab("coefficients")+
		ggplot2::ylab("") +
		mytheme
		print (pl)
	dev.off()
}
