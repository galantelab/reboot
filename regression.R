#!/usr/bin/env Rscript

#Time counter
start_time <- Sys.time()

######Argument parser######
option_list <- list(
  optparse::make_option(c("-I", "--filein"), action = "store", type = 'character', dest = "fname",
                        help = "Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival parameters"),
  optparse::make_option(c("-O", "--outprefix"), action = "store", type = 'character', dest = "out", default = "reboot",
                        help = "Output file prefix. Default: reboot"),
  optparse::make_option(c("-B", "--bootstrap"), action = "store", type = "integer", dest = "boot_iter", default = 1,
                        help = "Number of iterations for bootstrap simulation (int). Default: 1"),
  optparse::make_option(c("-G", "--groupsize"), action = "store", type = 'integer', dest = "nel", default = 10,
                        help = "Number of genes/transcripts to be selected in each bootstrap simulation (int). Default: 10"),
  optparse::make_option(c("-P", "--pcentfilter"), action = "store", type = 'double', dest = "pf", default = 0.3,
                        help = "Percentage of correlated gene/transcript pairs allowed in each iteration (double). Default: 0.3"),
  optparse::make_option(c("-V", "--varfilter"), action = "store", type = 'double', dest = "var", default = 0.01,
                        help = "Minimum normalized variance (0-1) required for each gene/transcript among samples (double). Default: 0.01"),
  optparse::make_option(c("-T", "--ty"), action = "store", type = 'character', dest = "ty", default = "gene",
                        help = "Declare which type of transcriptome data to be analyzed: gene or transcript (character). Default: gene"),
  optparse::make_option(c("-F", "--force"), action = "store_true", dest = "force", default = FALSE,
                        help = "To force overcome follow up variance filter and/or proportion filter for survival status (<20%), choose -F (boolean). Default: FALSE"),
  optparse::make_option(c("-N", "--ncores"), action = "store", type = "integer", dest = "n_cores", default = 1,
                        help = "Number of CPU cores to use for parallel bootstrap (int). Default: 1"),
  optparse::make_option(c("-S", "--seed"), action = "store", type = "integer", dest = "seed", default = 123,
                        help = "Random seed for reproducibility (int). Default: 123")
)

#Change variables
opo <- optparse::OptionParser(option_list = option_list, add_help_option = T)
in_object <- optparse::parse_args(opo)
logname <- in_object$out
outname <- paste(in_object$out, "_signature.txt", sep = "")
outplot <- in_object$out
force <- in_object$force
ty <- in_object$ty
seed <- in_object$seed

 
######read input file######
full_data <- read.table(in_object$fname, header = T, row.names = 1, check.names = F)
colnames(full_data) <- gsub("-","__",colnames(full_data))
if (ty == "gene") {
	bar <- 0.0035
} else if (ty == "transcript") {
	bar <- 0.011
} else {
	cat("Type has to be either gene or transcript\n")
	q(status = 0)
}

####Setting plot theme#####
mytheme <- ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                 text = ggplot2::element_text(face = "plain", colour = "black"),
                 axis.text = ggplot2::element_text(face = "bold", colour = "black"),
                 legend.text = ggplot2::element_text(colour = "black", face = "plain"),
                 legend.title = ggplot2::element_text(colour = "black", face = "bold"),
                 axis.ticks = ggplot2::element_line(colour = "black"),
                 axis.line = ggplot2::element_line(colour = "black"))

#######Log file##########
sink(file = paste(logname, ".log", sep = ''))

# To be added in the beginning of the script since this is intended to be the first lines of the ".log" file
cat("\n\n============================================================")
cat(" Make Signature ")
cat("============================================================\n\n")

scriptname <- "Regression analysis"
cat("Chosen parameters: ")
cat(paste(commandArgs(trailingOnly = T), collapse = " "))
cat("\n\n")

#####Schoenfeld test######
ph_assumptions <- function(full_data) {
  cat("Performing schoenfeld test\n\n")
	filt <- vector()
	attributes <- colnames(full_data[3:dim(full_data)[2]])
	
	for (i in attributes)
	{
	  phmodel <- survival::coxph(formula = formula(paste('survival::Surv(OS.time, OS)~', i)) , data = full_data)
	  if (!is.na(phmodel$coef))
	  {
	    tryCatch(
	      {
	        phmodel <- survival::coxph(formula = formula(paste('survival::Surv(OS.time, OS)~', i)) , data = full_data)
	        schoen <- survival::cox.zph(phmodel)
	        pval <- schoen$table[1,3]
	        if (pval > 0.05) {filt <- c(filt, i)}
	      }, warning = function(w){}, error = function(e){})
	  }
	}
	
	losers <- setdiff(attributes, filt)
	cat(length(losers)," columns not allowed by schoenfeld test: ",losers, "\n\n")
	
	return(full_data[,c(colnames(full_data)[1:2],filt)])
}

######Error check 1#######
numberfilter1 <- function(dataf, g, outname, outplot) {
	if ((ncol(dataf) - 2) < g)
	{
	  cat("The number of columns per group exceeds the number of columns", "\n", "\n")
	  cat("Performing single multivariate regression","\n","\n")
	  coemale <- regression(dataf, seed)
	  feature <- names(coemale)
	  coefficient <- unname(coemale)
	  
	  if (any(!(coefficient == 0)))
	  {
	    coemale <- cbind(feature,coefficient)
	    coemale <- gsub("__","-",coemale)
	    write.table(coemale, outname, sep = "\t", row.names = F, quote = F)
	    cat("Done", "\n")
	  } else {
	    cat("No signature found, all coefficients are equal 0", "\n")
	    sink()
	  }
	  q(status = 0)
	}
}

######Error check 2#######
numberfilter2 <- function(dataf, g, outname, outplot) {
  if (ncol(dataf) <= 2)
  {
    sink(file = paste(out, ".err", sep = ''), append = T)
    cat("No column was left after variance filter", "\n","\n")
    sink()
    q(status = 0)
  }
  
  if ((ncol(dataf) - 2) < g)
  {
    cat("The number of columns is lower than group size due to variance filter","\n","\n")
    cat("Performing single multivariate regression","\n","\n")
    coemale <- regression(dataf, seed)
    feature <- names(coemale)
    coefficient <- unname(coemale)
    
    if (any(!(coefficient == 0)))
    {
      coemale <- cbind(feature,coefficient)
      coemale <- gsub("__","-",coemale)
      write.table(coemale, outname, sep = "\t", row.names = F, quote = F)
      histogram(outplot,coemale)
      cat("Done", "\n")
    } else {
      cat("No signature found, all coefficients are equal 0", "\n")
    }
    q(status = 0)
  }
}

######Correlation filter######
corfun <- function(cmatrix, pf) {
  indexes <- c()
  pval <- c()
  names <- c()
  ngenes <- 3:(ncol(cmatrix) - 1)
  
  for (t in ngenes)
  {
    for (u in ((t + 1):(ncol(cmatrix))))
    {
      aux <- suppressWarnings(cor.test(x = cmatrix[,t], y = cmatrix[,u], method = 'spearman'))
      indexes <- c(indexes, round(aux$estimate, 3))
      pval <- c(pval, round(aux$p.value, 3))
      names <- c(names, paste(colnames(cmatrix)[t], colnames(cmatrix)[u], sep = "_"))
    }
  }
  
  if (((sum((indexes > 0.80) & (pval < 0.05))) / length(pval)) >= pf)
  {
    switch <- 1
    cat("This iteration was avoided due to correlation among columns. Sperman correlation values and p-values are respectively:", "\n")
    cat(paste(names, indexes, pval, sep = ":"), "\n", "\n")
  } else {
    switch <- 0
  }
  
  return(switch)
}

######Bootfunction######
bootstrapfun <- function(full_data, boot_iter, nel, outname, outplot, pf, bar, n_cores = 1) {
	cat("Starting bootstrap ", boot_iter, "iterations\n\n")
	k <- colnames(full_data[3:length(colnames(full_data))])
	v <- vector("list", length(k))
	yield <- hash::hash(k, v)

	collect_coefficients <- function(results_list) {
	  for (coemale in results_list) {
	    for (j in 1:length(coemale)) {
	      name <- names(coemale[j])
	      val  <- coemale[j][[1]]
	      eval(parse(text = paste("yield$", name, "<- c(yield$", name, ",", val, ")", sep = "")))
	    }
	  }
	}

	if (n_cores == 1) {
	  i <- 1
	  while (i <= boot_iter) {
	    cat("processing iteration: ", i, "\n\n")
	    if (ncol(full_data) == 3) {
	      cmatrix <- full_data
	    } else {
	      cmatrix <- subsample(full_data, nel, i)
	      if (nel > 1) {if (corfun(cmatrix, pf) == 1) {next}}
	    }
	    coemale <- regcall(cmatrix, nel, full_data)
	    if (!is.null(coemale)) {
	      i <- i + 1
	      for (j in 1:length(coemale)) {
	        name <- names(coemale[j])
	        val  <- coemale[j][[1]]
	        eval(parse(text = paste("yield$", name, "<- c(yield$", name, ",", val, ")", sep = "")))
	      }
	    }
	  }
	} else {
	  cat("Using ", n_cores, " cores for parallel bootstrap\n\n")

	  # One attempt per seed: subsample -> correlation check -> regression
	  run_attempt <- function(seed_val) {
	    if (ncol(full_data) == 3) {
	      cmatrix <- full_data
	    } else {
	      set.seed(seed_val)
	      shuffle  <- sample(colnames(full_data[, 3:ncol(full_data)]), size = nel, replace = FALSE)
	      cmatrix  <- cbind(full_data[, 1:2], subset(full_data, select = shuffle))
	      if (nel > 1 && corfun(cmatrix, pf) == 1) return(NULL)
	    }
	    tryCatch(
	      R.utils::withTimeout(regression(cmatrix, seed = seed), timeout = 60, onTimeout = "warning"),
	      warning = function(w) NULL,
	      error   = function(e) NULL
	    )
	  }

	  cl <- parallel::makeCluster(n_cores)
	  on.exit(parallel::stopCluster(cl), add = TRUE)
	  parallel::clusterExport(cl, varlist = c("full_data", "nel", "pf", "seed", "corfun", "regression"),
	                          envir = environment())
	  parallel::clusterEvalQ(cl, {
	    library(survival)
	    library(penalized)
	    library(R.utils)
	  })

	  valid_results <- list()
	  seed_counter  <- 1
	  cat("Bootstrap progress: 0/", boot_iter, "\n", sep = "")

	  while (length(valid_results) < boot_iter) {
	    remaining   <- boot_iter - length(valid_results)
	    batch_size  <- max(remaining, n_cores)
	    batch_seeds <- seed_counter:(seed_counter + batch_size - 1)
	    seed_counter <- seed_counter + batch_size

	    batch       <- parallel::parLapply(cl, batch_seeds, run_attempt)
	    valid_batch <- Filter(Negate(is.null), batch)
	    valid_results <- c(valid_results, valid_batch)
	    cat("Bootstrap progress: ", min(length(valid_results), boot_iter), "/", boot_iter, "\n", sep = "")
	  }

	  collect_coefficients(valid_results[1:boot_iter])
	  cat("\n")
	}

	#Processing result
	aux <- c(NULL, NULL, NULL)
	for (feature in hash::keys(yield)) {
	  coefficient <- suppressWarnings(eval(parse(text = paste("mean(yield$", feature, ")", sep = ""))))
	  sd          <- suppressWarnings(eval(parse(text = paste("sd(yield$", feature, ")", sep = ""))))
	  aux <- rbind(aux, cbind(feature, coefficient, sd))
	}

	tt <- as.data.frame(aux)
	if (any(!complete.cases(tt$coefficient))) {cat("NA coefficient found, increase coverage for a proper analysis", "\n")}
	tt <- dplyr::filter(tt, abs(as.numeric(as.character(coefficient))) >= bar)

	if (any(!(tt$coefficient == 0)) & dim(tt)[1] != 0) {
	  tt$feature <- gsub("__", "-", tt$feature)
	  write.table(tt, outname, sep = "\t", row.names = F, quote = F)
	  histogram(outplot, tt)
	  lolli(outplot, tt)
	} else {cat("No signature found, all coefficients are not significant", "\n")}
	cat("Done\n")
}

######Variance filter######
varfun <- function(cmatrix, var, file, force, out) {
	if (any(apply(cmatrix, 2, function(x) any(is.na(x)))))
	{
	  cat("Cheking NAs\n")
	  impu <- mice::mice(cmatrix, print = F)
	  cmatrix <- mice::complete(impu)
	}
	
  colnames(cmatrix)[1:2] <- c("OS", "OS.time")

	if (ncol(cmatrix) == 3)
	{
	  divisor <- max(cmatrix[3])
	  dividend <- cmatrix[3]
	  normalized <- divisor/dividend
	  variance <- var(normalized)
	  
	  if (variance < var)
	  {
	    cat("All columns rejected by variance filter", "\n", "\n")
	    q(status = 0)
	  }
	  cat("No columns rejected by variance filter", "\n", "\n")
	} else {
	  maxes <- matrix(apply(cmatrix[,3:ncol(cmatrix)], 2, max), nrow = 1)
	  
	  if (0 %in% maxes)
	  {
	    cat("Columns with only 0s found in ", file, ". Remove such columns and try again.", "\n")
	    q(status = 0)
	  }
	  
	  cat("Calculating normalized variances", "\n\n")
	  dividend <- dplyr::bind_rows(replicate(nrow(cmatrix), as.data.frame(maxes), simplify = F))
	  divisor <- cmatrix[,3:ncol(cmatrix)]
	  normalized <- divisor/dividend
	  variances <- apply(normalized,2,var)
	  filtered <- c()
	  losers <- c()
	  
	  for (i in (1:length(variances))) {if (variances[i] > var) {filtered <- c(filtered, i)}}
	  
	  if (class(filtered) == "NULL")
	  {
	    cat("All columns rejected by variance filter","\n","\n")
	    q(status = 0)
	  } else if (length(filtered) == length(variances)) {
	    cat("No columns rejected by variance filter", "\n", "\n")
	  } else {
	    losers <- names(variances[-filtered])
	    filtered <- filtered + 2
	    cmatrix <- cmatrix[, c(1,2,filtered)]
	    cat(length(losers)," columns with variance lower than ", var, " was removed from analysis: ",losers, "\n","\n")
	  }
	  
	  #Dealing with SO and SO time
	  if (!force)
	  {
	    OSstatus <- cmatrix[,1]
	    percentage <- sum(OSstatus)/length(OSstatus)
	    
	    if (percentage < 0.2 | percentage > 0.8)
	    {
	      sink(file = paste(out, ".err", sep = ''), append = T)
	      cat("Survival status proportion:", percentage, " is probably not enough to the analysis. \nDeath or recidive are alternative options for the analysis", "\n\n")	
	      cat("If you want to continue anyway, choose the flag F. \n\n")
	      sink()
	      q(status = 0)
	    }
	    
	    followup <- cmatrix[,2]
	    uplimit <- max(followup)	
	    normalized <- followup/uplimit
	    fvar <- var(normalized)
	    
	    if (fvar < var)
	    {
	      sink(file = paste(out, ".err", sep = ''), append = T)
	      cat("Follow up variance: ", var, " has not passed the variance test. \n\n")
	      cat("If you want to continue anyway, choose the flag F, or lower variance filter. \n\n")
	      sink()
	      q(status = 0)
	    }
	  }
	}
	
  return(cmatrix)
}

######Subsampling procedure#########
subsample <- function(full_data, nel, seed=i) {
  set.seed(seed)
  shuffle <- sample(colnames(full_data[,3:ncol(full_data)]), size = nel, replace = F)
  cat("Picked columns: ", shuffle, "\n", "\n")
  cmatrix <- cbind(full_data[,1:2], subset(full_data, select = shuffle))
}

######Regression time keeper#######
regcall <- function(cmatrix, nel, full_data) {
  tryCatch(R.utils::withTimeout(coemale <- regression(cmatrix, seed), timeout = 60, onTimeout = "warning"),
           warning = function(warning_condition) {
             cat("Regression time exceeded, you may consider changing variance and/or correlation filters. Trying again \n");
             cmatrix <- subsample(full_data, nel);
             regcall(cmatrix, nel, full_data)
           }, error = function(e){}
  )
}

######Regression procedure############
regression <- function(cmatrix, seed = seed) {
  set.seed(seed)	
  formula <- survival::Surv(OS.time, OS) ~ .
  fit1 <- penalized::profL1(formula, data = cmatrix, fold = 10, plot = F, trace = F)
  opt1 <- penalized::optL1(formula, data = cmatrix, fold = fit1$fold, trace = F)
  fit <- penalized::penalized(formula, data = cmatrix, lambda1 = opt1$lambda, trace = F)
  coemale <- penalized::coefficients(fit, "all")
  return(coemale)
}

#######lollipop plot#################
lolli <- function(out, tt) {
  cat("building ranking plot\n")
  fname <- paste(out, "_lollipop.pdf", sep = "")
  tt <- tt[complete.cases(tt), ]
  tt$coefficient <- as.numeric(as.character(tt$coefficient))
  tt <- dplyr::filter(tt,coefficient != 0)
  tt$feature = factor(tt$feature, levels = tt[order(tt$coefficient), "feature"])
  filter <- as.character(dplyr::top_n(tt, 10, abs(coefficient))$feature)
  tt <- dplyr::filter(tt, feature %in% filter)
  tt <- dplyr::mutate(tt, name = forcats::fct_reorder(feature, coefficient))
  
  pdf(fname)					
  ll <- ggplot2::ggplot(tt, ggplot2::aes(x = reorder(feature, coefficient), y = coefficient)) +
    ggplot2::geom_segment(ggplot2::aes(x = feature, xend = feature, y = 0, yend = coefficient), linewidth = 2, color = "grey") +
    ggplot2::geom_point(size = 4, colour = "#1c9099") + mytheme +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  print(ll, newpage = F)
  dev.off()
}

#######Histogram plot#################
histogram <- function(out, tt) {
  cat("building histogram plot\n")
  fname <- paste(out, "_hist.pdf", sep = "")
  tt <- tt[complete.cases(tt), ]
  tt$coefficient <- as.numeric(as.character(tt$coefficient))
  tt <- dplyr::filter(tt, coefficient != 0)
  
  pdf(fname)
  pl <- ggplot2::ggplot(tt, ggplot2::aes(x = coefficient)) +
    ggplot2::geom_histogram(binwidth = 0.005, alpha = 1, position = "identity") +
    ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::xlab("coefficients") + ggplot2::ylab("") + mytheme
  print(pl)
  dev.off()
}

####Main####

#Perform variance filter
full_data <- varfun(full_data, in_object$var, in_object$fname, force, logname)

#Perform schoenfeld tests
full_data <- ph_assumptions(full_data)

#check file error
numberfilter1(full_data, in_object$nel, outname, outplot)

#check file error 2
numberfilter2(full_data, in_object$nel, outname, outplot)

#Perform Regression
bootstrapfun(full_data, in_object$boot_iter, in_object$nel, outname, outplot, in_object$pf, bar, in_object$n_cores)

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
