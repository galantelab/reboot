start_time <- Sys.time()

suppressMessages(library("optparse"))

option_list <- list(

make_option(c("-I", "--filein"), action="store",
type='character', dest = "exp_file", help="Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival paramenters"),

make_option(c("-O", "--outprefix"), action="store",
type='character', dest = "out", default = "reboot", help="Output file prefix. Default: reboot"),

make_option(c("-S", "--signature"), action="store",
type='character', dest = "sig_file", help="Tab separated values (tsv) file containing a set of genes/transcripts and corresponding cox coefficients"),

make_option(c("-M", "--multivariate"), action="store",
type='logical', dest = "type", default = FALSE, help="If clinical variables should be included, choose -M. This option is tied with -C option"),

make_option(c("-C", "--clinical"), action="store",
type='character', dest = "clin_file", default = "", help="Tab separated values (tsv) file containing binary categorical variables only. Required if -M option is chosen"),

make_option(c("-R", "--roc"), action="store",
type='logical', dest = "roc_curve", default = FALSE, help="If If genetic score should be categorized according to a ROC curve instead of median, choose -R"))

opo <- OptionParser(option_list=option_list, add_help_option = T)
in_object <- parse_args(opo)

#Check if clinical data is provided in case type == TRUE
if(in_object$type){
  if(in_object$clin_file == ""){
    cat("Insert clinical variables file\n")
    q(status=0)
  }
}

#Change variables
exp_file = in_object$exp_file
out = in_object$out
sig_file = in_object$sig_file
type = in_object$type
clin_file = in_object$clin_file
roc_curve = in_object$roc_curve

####### Log file ##########
sink(file = paste(out, ".log", sep=''), append=T)

#Load libraries
suppressMessages(library("survcomp"))
suppressMessages(library("survival"))
suppressMessages(library("survminer"))
suppressMessages(library("OptimalCutpoints"))
suppressMessages(library("survivalROC"))
suppressMessages(library("forestmodel"))
suppressMessages(library("sjstats"))
suppressMessages(library("data.table"))
suppressMessages(library("plyr"))
suppressMessages(library("dplyr"))

cat("\n\n============================================================")
cat(" Apply Signature ")
cat("============================================================\n\n")

cat("Chosen parameters: ")
cat(paste(commandArgs(trailingOnly = T), collapse = " "))
cat("\n\n")

#Print log message
cat("Checking provided files...\n")

#Check if provided files exist
if(!file.exists(exp_file)){
  cat(paste("Error: File ", exp_file, " does not exist. Please provide a valid file.\n",sep=""))
  quit(save="no")
}
if(!file.exists(sig_file)){
  cat(paste("Error: File ", sig_file, " does not exist. Please provide a valid file.\n",sep=""))
  quit(save="no")
}
if(type & clin_file!="" & !file.exists(clin_file)){
  cat(paste("Error: File ", clin_file, " does not exist. Please provide a valid file.\n",sep=""))
  quit(save="no")
}

#Read provided files
signature = read.table(sig_file, header=T, check.names=F, stringsAsFactors=F)
data = read.table(exp_file, header=T, check.names=F, stringsAsFactors=F)

if(type & clin_file!=""){
  clin = read.delim(clin_file, header=T, row.names=1)
}

#Check number of columns in signature file
if(ncol(signature)!=2){
  cat(paste("Error: ", sig_file, " should have 2 columns. Please check the manual for more information.\n",sep=""))
  quit(save="no")
}

#Check signature file format
if(!is.character(signature[,1])){
  cat(paste("Error: Invalid format for column 1 of file ", sig_file, ". It should be a list of features (characters).\n",sep=""))
  quit(save="no")
}

#Check signature format
if(!is.numeric(signature[,2])){
  cat(paste("Error: Invalid format for column 2 of file ", sig_file, ". It should be a list of Cox coefficients (numeric values).\n",sep=""))
  quit(save="no")
}

colnames(signature) = c("feature", "coefficient")
colnames(data)[1:3] = c("sample", "OS", "OS.time")

#Replace column names
signature$feature = gsub("-","_",signature$feature)
colnames(data) = gsub("-","_",colnames(data))

#Check if there is at least one valid feature in signature
signature = signature[!is.na(signature$coefficient) & signature$coefficient!=0,]
if(nrow(signature)==0){
  cat(paste("Error: No valid features to be tested in the signature.\n",sep=""))
  quit(save="no")
}

#Check if all genes in signature are present on data file
data = data[,grep(paste("^sample$","^OS$","^OS.time$",paste("^",signature$feature,"$",collapse="|",sep=""),sep="|"),colnames(data))]
if(ncol(data) != (nrow(signature)+3)){
  if(ncol(data) > (nrow(signature)+3)){
    cat(paste("Error: Duplicated features are not allowed. Please check the provided files.\n",sep=""))
    quit(save="no")
  }
  if(ncol(data) < (nrow(signature)+3)){
    cat(paste("Error: Not all features are present in both files. Please check the provided files.\n",sep=""))
    quit(save="no")
  }
}

#Reorder colnames of data according to signature order
target = c("sample","OS","OS.time",signature$feature)
data = data[,match(target, colnames(data))]

#Create rownames
rownames(data) = data$sample
data = data[,2:ncol(data)]

#Create signature score for each sample
tmp = as.matrix(data[,3:ncol(data)])
if(nrow(signature)>1){ tmp = as.data.frame(tmp%*%diag(signature$coefficient)) }else{ tmp = as.data.frame(tmp*signature$coefficient) }
tmp$score <- apply(tmp, 1, function(x) sum(x))
tmp = tmp[,ncol(tmp),drop=F]
data = cbind(data[,1:2],tmp)

#Edited "ggcoxzph" function to add the global Schoenfeld Test p-value
reboot_ggcoxzph <- function(fit, resid = T, se = T, df = 4, nsmo = 40, var, point.col = "red", point.size = 1,
                            point.shape = 19, point.alpha = 1, caption = NULL, ggtheme = theme_survminer(), ...){
  x <- fit
  if(!methods::is(x, "cox.zph"))
    stop("Can't handle an object of class ", class(x))
  
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
  if (qmat$rank < df)
    stop("Spline fit is singular, try a smaller degrees of freedom")
  if (se) {
    bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
    xtx <- bk %*% t(bk)
    seval <- d * ((pmat %*% xtx) * pmat) %*% rep(1, df)
  }
  ylab <- paste("Beta(t) for", dimnames(yy)[[2]])
  if (missing(var))
    var <- 1:nvar
  else {
    if (is.character(var))
      var <- match(var, dimnames(yy)[[2]])
    if (any(is.na(var)) || max(var) > nvar || min(var) <
        1)
      stop("Invalid variable requested")
  }
  if (x$transform == "log") {
    xx <- exp(xx)
    pred.x <- exp(pred.x)
  }
  else if (x$transform != "identity") {
    xtime <- as.numeric(dimnames(yy)[[1]])
    indx <- !duplicated(xx)
    apr1 <- approx(xx[indx], xtime[indx], seq(min(xx), max(xx),
                                              length = 17)[2 * (1:8)])
    temp <- signif(apr1$y, 2)
    apr2 <- approx(xtime[indx], xx[indx], temp)
    xaxisval <- apr2$y
    xaxislab <- rep("", 8)
    for (i in 1:8) xaxislab[i] <- format(temp[i])
  }
  plots <- list()
  lapply(var, function(i) {
    invisible(round(x$table[i, 3],4) -> pval)
    invisible(round(x$table[nrow(x$table), 3],4) -> global)
    ggplot() + labs(title = paste0('Global Schoenfeld Test p: ', global),
                    subtitle = paste0('Individual Schoenfeld Test p: ', pval)) +
      ggtheme + theme(plot.title = element_text(hjust = .5, vjust = .5, face = "bold", margin = margin(0, 0, 10, 0)),
                      plot.subtitle = element_text(hjust = 0, vjust = .5, face = "plain", margin = margin(10, 0, 10, 0))) -> gplot
    y <- yy[, i]
    yhat <- as.vector(pmat %*% qr.coef(qmat, y))
    if (resid)
      yr <- range(yhat, y)
    else yr <- range(yhat)
    if (se) {
      temp <- as.vector(2 * sqrt(x$var[i, i] * seval))
      yup <- yhat + temp
      ylow <- yhat - temp
      yr <- range(yr, yup, ylow)
    }
    if (x$transform == "identity") {
      gplot + geom_line(aes(x=pred.x, y=yhat)) +
        xlab("Time") +
        ylab(ylab[i]) +
        ylim(yr) -> gplot
    } else if (x$transform == "log") {
      gplot + geom_line(aes(x=log(pred.x), y=yhat)) +
        xlab("Time") +
        ylab(ylab[i]) +
        ylim(yr)  -> gplot
    } else {
      gplot + geom_line(aes(x=pred.x, y=yhat)) +
        xlab("Time") +
        ylab(ylab[i]) +
        scale_x_continuous(breaks = xaxisval,
                           labels = xaxislab) +
        ylim(yr)-> gplot
    }
    
    if (resid)
      gplot <- gplot + geom_point(aes(x = xx, y =y),
                                  col = point.col, shape = point.shape, size = point.size, alpha = point.alpha)
    
    if (se) {
      gplot <- gplot + geom_line(aes(x=pred.x, y=yup), lty = "dashed") +
        geom_line(aes( x = pred.x, y = ylow), lty = "dashed")
    }
    
    ggpubr::ggpar(gplot, ...)
    
    
  }) -> plots
  names(plots) <- var
  class(plots) <- c("ggcoxzph", "ggsurv", "list")
  
  if("GLOBAL" %in% rownames(x$table)) # case of multivariate Cox
    global_p <- x$table["GLOBAL", 3]
  else global_p <- NULL # Univariate Cox
  attr(plots, "global_pval") <- global_p
  attr(plots, "caption") <- caption
  plots
}

#Function to test Cox Proportional Assumptions (Schoenfeld Test)
test_ph_assumptions <- function(model_object, covariates, is_multi)
{
  test.ph <- cox.zph(model_object)
  
  if(is_multi){

    tmp_covariates <- c()
    covariates = c("score", covariates)
    new_covariates = colnames(test.ph$y)
  
    #Match string to get from 'new covariates' the original 'covariates' name
    for(var in new_covariates)
    {
      for(var2 in covariates)
      {
        if(grepl(var2, var)){tmp_covariates <- append(tmp_covariates, var2)}
      }
    }
  } else{

    tmp_covariates <- c("score")

  }
  
  rownames(test.ph$table) <- c(tmp_covariates, "GLOBAL")
  colnames(test.ph$y) <- tmp_covariates
  
  if(is_multi) {
    phplot <- reboot_ggcoxzph(fit = test.ph)
  } else {
    phplot <- ggcoxzph(fit = test.ph)
  }
  
  pdf(file = paste(out, "_ph_assumptions_plot.pdf", sep=""))
  for (plot in phplot) {
    print(plot, newpage = T)
  }
  garbage = dev.off()
  
  tmp_df <- as.data.frame(test.ph$table)
  pvalue <- tmp_df[nrow(tmp_df),3]
  return(pvalue)

}

#Function to calculate ROC curve
cutoff_ROC <- function(dataset, auc_val, filename, plot)
{
  if (auc_val < .5)
  {
    optimal.cutpoint <- optimal.cutpoints(X = "score",
                                          status = "OS",
                                          tag.healthy = 0, methods = "Youden",
                                          data = dataset,
                                          control = control.cutpoints(),
                                          #ci.fit = TRUE,
                                          direction = ">")
  }
  else
  {
    optimal.cutpoint <- optimal.cutpoints(X = "score",
                                          status = "OS",
                                          tag.healthy = 0, methods = "Youden",
                                          data = dataset,
                                          #ci.fit = TRUE,
                                          control = control.cutpoints())
  }
  
  var_cutoff <- as.numeric(optimal.cutpoint$Youden$Global$optimal.cutoff$cutoff)
  
  if (length(var_cutoff) >= 2)
  {
    var_cutoff <- var_cutoff[1]
  }
  else
  {
    var_cutoff <- var_cutoff
  }
  
  if (plot){
    pdf(filename)
    roc_curve_plot <- plot.optimal.cutpoints(x = optimal.cutpoint, legend = T, which = c(1), col = "blue", bg = "white")
    garbage = dev.off()
  }
  
  return(var_cutoff)
}

#If some clinical file is provided
if(type & clin_file != ""){

  #Check if all samples in exp_file are present on clin_file
  clin = merge(data,clin,by=0)
  if(nrow(data) != nrow(clin)){
    if(nrow(data) > nrow(clin)){
      cat(paste("Error: Not all samples from file ", exp_file," are provided in file ", clin_file,". Please check the provided files.\n",sep=""))
      quit(save="no")
    }
    if(nrow(data) < nrow(clin)){
      cat(paste("Error: Not all samples from file ", clin_file," are provided in file ", exp_file,". Please check the provided files.\n",sep=""))
      quit(save="no")
    }
  }

  #Change rownames of clinical data
  rownames(clin) = clin$Row.names
  clin = clin[,2:ncol(clin)]

  #Print log message
  cat("Done\n")
  cat("\n")

  #Transform score to categorical with ROC curve or with median value
  if(roc_curve){
    #This value can be changed. By default, it uses the median followup time
    cutoff = median(x = clin$OS.time, na.rm = TRUE)
    nobs <- nrow(clin)
    
    roc = survivalROC(Stime = clin$OS.time, status = clin$OS, marker = clin$score, predict.time = cutoff, method = "NNE", span = 0.25*nobs^(-0.20))
    score_cutoff = cutoff_ROC(dataset = clin, auc_val = as.numeric(roc$AUC), plot = F)
    clin$score = ifelse(clin$score > score_cutoff, "high", "low")
  } else {
    clin$score = ifelse(clin$score > median(clin$score), "high", "low")
    }

  #Check if all columns have categorical variables with only 2 factors
  for(i in 3:ncol(clin)){

    if(nlevels(as.factor(clin[,i])) != 2){
      cat(paste("Error: Variable \"", colnames(clin)[i], "\" in file ", clin_file," does not have 2 categories. Please check the provided file.\n",sep=""))
      quit(save="no")
    }
  }
} else{

  #Print log message
  cat("Done\n\n")

}

#Print log message
cat("Creating signature score...\n")
cat("Done\n\n")

#Create function to run log-rank test for score signatures
logrank.test <- function(dat,filename){

  #Transform score to categorical with ROC curve or with median value
  if(roc_curve){
    #This value can be changed. By default, it uses the median followup time
    cutoff = median(x = dat$OS.time, na.rm = TRUE)
    nobs <- nrow(dat)

    cat("Generating ROC curve...\n")    
    roc = survivalROC(Stime = dat$OS.time, status = dat$OS, marker = dat$score, predict.time = cutoff, method = "NNE", span = 0.25*nobs^(-0.20))
    score_cutoff = cutoff_ROC(dataset = dat, auc_val = as.numeric(roc$AUC), filename = paste(out, "_ROC.pdf", sep=""), plot = T)
    dat$score = ifelse(dat$score > score_cutoff, "high", "low")
    cat("Done\n\n")

  } else {
    dat$score = ifelse(dat$score > median(dat$score), "high", "low")
  }

  #Test proportional hazards assumptions
  cat("Testing proportional hazards assumption (signature score)...\n")
  uni_model = coxph(formula = formula(paste('Surv(OS.time, OS) ~ score')) , data = dat)
  checkPH <- test_ph_assumptions(model_object = uni_model, covariates = "NULL", is_multi = F)
  if (checkPH <= .05) {
    cat(paste("Warning: Proportional Hazards Assumptions not met (p = ", round(x = checkPH, digits = 4), "). Check plot: ",
              paste("'", out, "_ph_assumptions_plot.pdf'\n", sep=""), sep = ""))
  } else {
    cat(paste("Proportional Hazards Assumptions met (p = ", round(x = checkPH, digits = 4), ").\n", sep = ""))
  }
  cat("Done\n\n")

  cat("Running log-rank test for signature score...\n")

  tryCatch({

    #Run log rank test
    hr_model = hazard.ratio(x = dat$score, surv.time = dat$OS.time, surv.event = dat$OS, alpha = .05, method.test = "logrank", na.rm = T)

    #Extract result fields: p-value, hazard ratio(95% confidence interval) and coefficient
    hazard.ratio = paste(round(hr_model$hazard.ratio,4), " (95% CI, ", round(hr_model$lower,4), " - ", round(hr_model$upper,4), ")", sep="")
    coef = hr_model$coef
    log.rank.pvalue = hr_model$p.value

    #Get number of high/low expression samples
    low.high.samples = paste(nrow(dat[dat$score=="low",]),nrow(dat[dat$score=="high",]),sep="/")

    #Run survfit to get median survival of groups and 95% confidence interval
    fit = survfit(Surv(OS.time, OS) ~ score, data=dat)
    median.survival.low = paste(summary(fit)$table[2,7], " (95% CI, ", summary(fit)$table[2,8]," - ", summary(fit)$table[2,9], ")" ,sep="")
    median.survival.high = paste(summary(fit)$table[1,7], " (95% CI, ", summary(fit)$table[1,8]," - ", summary(fit)$table[1,9], ")" ,sep="")

    #Add information about prognosis
    if(!is.na(log.rank.pvalue) & log.rank.pvalue<0.05 & coef<0){
      prognosis = "better"
    } else if(!is.na(log.rank.pvalue) & log.rank.pvalue<0.05 & coef>0){
      prognosis = "worse"
    } else{
      prognosis = "----"
    }

    #Create dataframe with result to be used in multivariate analysis
    result = data.frame(feature="score", coefficient=round(coef,4), hazard.ratio=hazard.ratio, log.rank.pvalue=round(log.rank.pvalue,4),
                        low.high.samples=low.high.samples, median.survival.low=median.survival.low, median.survival.high=median.survival.high, prognosis=prognosis)

    #Write log-rank result to file
    write.table(result, filename, row.names=F, col.names=T, quote=F, sep="\t")

    return(result)

  }, warning = function(w){}, error = function(e){})

}

#Run log-rank test for score signatures
res_logrank = logrank.test(data, paste(out, "_logrank.txt",sep=""))

#Print log messages
cat("Done\n\n")
cat("Generating Kaplan-Meier curve...\n")

#Create function to plot survival curve
logrank.plot <- function(dat,filename){

  #Get cutoff (either median or ROC) score of data signature AND transform data to categorical
  if(roc_curve){
    #This value can be changed. By default, it uses the median followup time
    cutoff = median(x = dat$OS.time, na.rm = TRUE)
    nobs <- nrow(dat)
    
    roc = survivalROC(Stime = dat$OS.time, status = dat$OS, marker = dat$score, predict.time = cutoff, method = "NNE", span = 0.25*nobs^(-0.20))
    score_cutoff = cutoff_ROC(dataset = dat, auc_val = as.numeric(roc$AUC), plot = F)
    sig_value = round(score_cutoff,2)
    dat$score = ifelse(dat$score > score_cutoff, "high", "low")
  } else {
    sig_value = round(median(dat$score),2)
    dat$score = ifelse(dat$score > median(dat$score), "high", "low")
  }

  #If there is at least two groups to be compared, create Kaplan-meier curve
  if(nlevels(as.factor(dat$score))>1){

    pdf(filename)
    pp = ggsurvplot(survfit(Surv(OS.time, OS) ~ score, data=dat), risk.table=TRUE, tables.theme = theme_cleantable(), tables.y.text = F, tables.height=0.2, pval=paste("p =",format(res_logrank[,4],scientific=T),sep=" "),
                    font.legend=16, font.x=22, font.y=22, font.tickslab=8, pval.size=6, pval.coord=c(0,0.05), title= "", legend = c(0.7, 0.9), legend.title="", break.time.by=200,
                    censor=T, legend.labs = c(paste("score>",sig_value, sep=""), paste("score<=",sig_value, sep="")), data=dat) + xlab("Survival time (days)")
    print(pp, newpage = FALSE)
    garbage = dev.off()

  } else{
    cat(paste("Warning: Data could not be partitioned into low/high scores.\n",sep=""))
  }
}

#Make survival plot (log-rank test for score signatures)
logrank.plot(data, paste(out, "_km_plot.pdf",sep=""))

#Print log messages
cat("Done\n\n")

#Create function to run univariate Cox-regression for each provided clinical parameter
univCox.test <- function(dat, covariates){

  univ_formulas = sapply(covariates, function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
  univ_models = lapply( univ_formulas, function(x){coxph(x, data = dat)})
  univ_results <- lapply(univ_models, function(x){ 

                       coef = as.data.frame(coef(summary(x)))
                       ci = as.data.frame(summary(x)$conf.int)
                       res = merge(coef,ci,by=0)
                       res = res[,c(1:3,6,9:10)]
                       colnames(res) = c("variable", "coefficient", "hr", "Cox.pvalue", "lower.ci", "upper.ci")
                       res$hazard.ratio = paste(round(res$hr,4), " (95% CI, ", round(res$lower.ci,4), " - ", round(res$upper.ci,4), ")", sep="")
                       res$Cox.pvalue = round(res$Cox.pvalue,4)
                       res$prognosis = "----"
                       if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.2 & res$coefficient<0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.2 & !is.na(res$coefficient) & res$coefficient<0,]$prognosis <- "better"}
                       if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.2 & res$coefficient>0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.2 & !is.na(res$coefficient) & res$coefficient>0,]$prognosis <- "worse"}
                       res = res[,c(1,7,4,8)]
                       return(res)
                })
  univ_res = do.call("rbind", univ_results)

  return(univ_res)
}

#Bootstrap Resampling
my_bootstrap_method <- function(raw_df, boot_df, boot_sample)
{
  boot_sample <- boot_sample + 1
  boot_vec <- boot_df$strap[[boot_sample]]$id
  bootstrap_df <- as.data.frame(raw_df[boot_vec,])
  
  #Check if all columns have categorical variables with only 2 factors
  for(i in 3:ncol(bootstrap_df))
  {
    if(nlevels(as.factor(bootstrap_df[,i])) != 2)
    {
      bootstrap_df <- "repeat"
      break
    } else {
      myLevels <- levels(as.factor(bootstrap_df[,i]))
      level1 <- (nrow(bootstrap_df[bootstrap_df[[i]] == myLevels[[1]],]) / nrow(bootstrap_df))
      if (level1 < .2 | level1 > .8)
      {
        bootstrap_df <- "repeat"
        break
      }
    }
  }
  
  return(bootstrap_df)
}

#Make forest plot for multivariate analysis
generate_forest_plot <- function(model_object, filename)
{
  pdf(filename)
  forest_plot <- suppressWarnings(forest_model(model = model_object, exponentiate = T, factor_separate_line = F, recalculate_width = T, recalculate_height = T))
  suppressWarnings(print(forest_plot, newpage = FALSE))
  garbage = dev.off()
}

#Merge tmp tables generated in multivariate cox regression with bootstrap
merge_tables <- function(tables_list, covariates)
{
  cox_tables <- list()
  iterator <- 0
    
  for (i in 1:length(tables_list))
  {
    #Remove rows (variables) if p-value > .05
    table <- tables_list[[i]]
    table <- table[!is.na(table$prognosis),]
    iterator <- iterator + 1
    tmp_df <- count(table, variable)
    cox_tables[[iterator]] <- tmp_df
  }
  
  tmp_table <- rbindlist(l = cox_tables, use.names = T, fill = T, idcol = "Unique ID")
  colnames(tmp_table) <- c("Unique ID", "Co-Variables", "Frequency")
  
  final_table <- count(tmp_table, `Co-Variables`)
  colnames(final_table) <- c("Co-Variables", "Frequency")
  
  tmp_covariates <- c()
  covariates = c("score", covariates)
  new_covariates = final_table[[1]]
  
  #Match string to get from 'new covariates' the original 'covariates' name
  for(var in new_covariates)
  {
    for(var2 in covariates)
    {
      if(grepl(var2, var)){tmp_covariates <- append(tmp_covariates, var2)}
    }
  }
  
  final_table[[1]] = tmp_covariates
  
  return(final_table)
}

#Make histograms of frequency of variables in each bootstrap multivariate model
barplot_co_variables <- function(plot_df, filename, covariates)
{
  tmp_covariates <- c()
  covariates = c("score", covariates)
  new_covariates = plot_df[[1]]

  #Match string to get from 'new covariates' the original 'covariates' name
  for(var in new_covariates)
  {
    for(var2 in covariates)
    {
      if(grepl(var2, var)){tmp_covariates <- append(tmp_covariates, var2)}
    }
  }

  plot_df$var = tmp_covariates

  pdf(filename)
  final_plot <- ggplot(data = plot_df, aes(x = reorder(plot_df[[3]], -plot_df[[2]]), y = plot_df[[2]])) +
    geom_bar(stat = "identity", width = 0.5, color = "black", fill = "black") +
    xlab("") + ylab("Frequency (%)") +
    geom_segment(aes(x = .5, y = 25, xend = (nrow(plot_df) + .5), yend = 25), color = "red", linetype = "dashed", size = .5) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.text = element_text(face = "plain", colour = "black"), 
    legend.text = element_text(colour = "black", face = "plain"), axis.ticks = element_line(colour = "black"), axis.line = element_line(colour = "black"), 
    panel.grid.major.x = element_blank(), text = element_text(face = "plain", colour = "black")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20))
  print(final_plot, newpage = FALSE)
  garbage = dev.off()

}

#Create function to run multivariate Cox regression for score + clinical parameters
multiCox.test <- function(dat, univ_result, logrank_result, uni_covariates){
  
  roc_curve2 <<- F
  filters <<- T

  #Select only relevant clinical parameters for multivariable Cox regression (p<0.2)
  covariates = c("score",rownames(univ_result[(!is.na(univ_result$Cox.pvalue) & univ_result$Cox.pvalue<0.2),]))
  
  tryCatch({

    if (roc_curve){
      
      #Remove row if "NA" in at least one column (variable)
      counter <- ncol(dat)
      
      tmp_table <- nrow(dat[complete.cases(dat),])
      sample_cutoff <- .7 * nrow(dat)
      temporary_table <- dat
      backup_dat <- dat
      
      # Ensure table without 'NAs' has at least 70% of the samples from the original dataframe
      while (tmp_table < sample_cutoff)
      {
        
        if (counter < 5){
          cat("\tWarning: Minimum number of co-variables (3) not achieved. Please check variables provided in clinical file. ")
          cat("Performing multivariate regression without bootstrap resampling...\n\n")
          break
        } else {
          
          remove_var_vec <- c()
          remove_index <- 0
          
          for (col in names(temporary_table))
          {
            tmp_var_size <- as.vector(is.na(temporary_table[[col]]))
            tmp_var_size2 <- length(tmp_var_size[tmp_var_size == T])
            
            remove_index <- remove_index + 1
            remove_var_vec[remove_index] <- tmp_var_size2
          }
          
          remove_col <- which.max(remove_var_vec)
          
          if (remove_col == length(remove_var_vec))
          {
            temporary_table <- temporary_table[,c(1:(remove_col-1))]
          } else {
            temporary_table <- temporary_table[,c(1:(remove_col-1),(remove_col+1):ncol(temporary_table))]
          }
          
          tmp_table <- nrow(temporary_table[complete.cases(temporary_table),])
          dat <- temporary_table[complete.cases(temporary_table),]
          counter <- counter - 1
        }
      }
      
      n_cols_bef <- ncol(dat)
      tmp_cols <- c()

      # Ensure less abundant category of each variable in table is at least 20%
      for (col in names(dat))
      {
        var <- table(dat[[col]])
        var1 <- var[[1]]
        var2 <- var[[2]]

        if (var1 <= var2)
        {
          prop <- (var1 / (var1 + var2))

          if (prop < .2)
          {
            next
          } else
          {
            tmp_cols <- append(tmp_cols, col)
          }
        } else
        {
          prop <- (var2 / (var1 + var2))

          if (prop < .2)
          {
            next
          } else
          {
            tmp_cols <- append(tmp_cols, col)
          }
        }
      }
      
      dat <- dat[, tmp_cols]
      
      n_cols_aft <- ncol(dat)
      n_cols_diff <- n_cols_bef - n_cols_aft
      
      counter <- counter - n_cols_diff
      barPlot_counter <- counter
      
      if (counter < 5){
        
        filters <<- F
        
        cat("\tWarning: Minimum number of co-variables (3) not achieved. Please check variables provided in clinical file. ")
        cat("Performing multivariate regression without bootstrap resampling...\n\n")
        dat <- backup_dat
        
        #Run multivariate Cox
        res = coxph(formula = formula(paste('Surv(OS.time, OS)~', paste(covariates,collapse=" + "))) , data = dat)
        
        #Extract results
        coef = as.data.frame(coef(summary(res)))
        ci = as.data.frame(summary(res)$conf.int)
        res = merge(coef,ci,by=0)
        res = res[,c(1:3,6,9:10)]
        
        colnames(res) = c("variable", "coefficient", "hr", "Cox.pvalue", "lower.ci", "upper.ci")
        res$hazard.ratio = paste(round(res$hr,4), " (95% CI, ", round(res$lower.ci,4), " - ", round(res$upper.ci,4), ")", sep="")
        res$Cox.pvalue = round(res$Cox.pvalue,4)
        res$prognosis = NA
        if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & res$coefficient<0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & !is.na(res$coefficient) & res$coefficient<0,]$prognosis <- "better"}
        if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & res$coefficient>0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & !is.na(res$coefficient) & res$coefficient>0,]$prognosis <- "worse"}
        res = res[,c(1,7,4,8)]
      } else {
        
        # Update covariates
        update_covariates <- c()
        for (col in 3:ncol(dat))
        {
          tmp_colName <- colnames(dat[col])
          if (tmp_colName %in% covariates)
          {
            update_covariates <- append(x = update_covariates, values = tmp_colName)
          }
        }
        
        if (length(update_covariates) < 3){
          
          filters <<- F
          
          cat("\tWarning: Minimum number of co-variables (3) not achieved. Please check variables provided in clinical file. ")
          cat("Performing multivariate regression without bootstrap resampling...\n\n")
          dat <- backup_dat
          
          #Run multivariate Cox
          res = coxph(formula = formula(paste('Surv(OS.time, OS)~', paste(covariates,collapse=" + "))) , data = dat)
          
          #Extract results
          coef = as.data.frame(coef(summary(res)))
          ci = as.data.frame(summary(res)$conf.int)
          res = merge(coef,ci,by=0)
          res = res[,c(1:3,6,9:10)]
          
          colnames(res) = c("variable", "coefficient", "hr", "Cox.pvalue", "lower.ci", "upper.ci")
          res$hazard.ratio = paste(round(res$hr,4), " (95% CI, ", round(res$lower.ci,4), " - ", round(res$upper.ci,4), ")", sep="")
          res$Cox.pvalue = round(res$Cox.pvalue,4)
          res$prognosis = NA
          if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & res$coefficient<0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & !is.na(res$coefficient) & res$coefficient<0,]$prognosis <- "better"}
          if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & res$coefficient>0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & !is.na(res$coefficient) & res$coefficient>0,]$prognosis <- "worse"}
          res = res[,c(1,7,4,8)]
        } else {
          
          covariates <- update_covariates
          
          roc_curve2 <<- T
          
          tables_list <- list()
          boot_backup_dat <- dat
          
          boot_tmp <- bootstrap(data = dat, n = 100, size = .6)
          boot_loop <- 0
          boot_control <- 0
          cat("\tBootstrap progress:\n\t")
          
          while(boot_loop < 100)
          {
            boot_loop <- boot_loop + 1
            dat <- my_bootstrap_method(raw_df = boot_backup_dat, boot_df = boot_tmp, boot_sample = boot_control)
            
            if (class(dat) != "data.frame") {
              boot_loop <- boot_loop - 1
              boot_control <- 0
              boot_tmp <- bootstrap(data = boot_backup_dat, n = (100 - boot_loop), size = .6)
            } else {
              boot_control <- boot_control + 1
              
              if (boot_loop == 100) {cat("# 100%\n\n")} else {cat("#")}
              
              #Run multivariate Cox
              res = coxph(formula = formula(paste('Surv(OS.time, OS)~', paste(covariates, collapse = " + "))) , data = dat)
              
              #Extract results
              coef = as.data.frame(coef(summary(res)))
              ci = as.data.frame(summary(res)$conf.int)
              res = merge(coef,ci,by=0)
              res = res[,c(1:3,6,9:10)]
              
              colnames(res) = c("variable", "coefficient", "hr", "Cox.pvalue", "lower.ci", "upper.ci")
              res$hazard.ratio = paste(round(res$hr,4), " (95% CI, ", round(res$lower.ci,4), " - ", round(res$upper.ci,4), ")", sep="")
              res$Cox.pvalue = round(res$Cox.pvalue,4)
              res$prognosis = NA
              if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & res$coefficient<0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & !is.na(res$coefficient) & res$coefficient<0,]$prognosis <- "better"}
              if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & res$coefficient>0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & !is.na(res$coefficient) & res$coefficient>0,]$prognosis <- "worse"}
              res = res[,c(1,7,4,8)]
              
              tables_list[[boot_loop]] <- res
            }
          }
          
          merged_table <<- merge_tables(tables_list = tables_list, covariates = uni_covariates)
          
          #Select only frequent (at least 50%) parameters for multivariable cox regression
          tmp_merged_table <- merged_table[merged_table$Frequency >= 25,]
          new_covariates <- as.vector(tmp_merged_table[[1]])
          tmp_covariates <- c()
          
          #Match string to get from 'new covariates' the original 'covariates' name
          for(var in new_covariates)
          {
            for(var2 in covariates)
            {
              if(grepl(var2, var)){tmp_covariates <- append(tmp_covariates, var2)}
            }
          }
          
          #Run multivariate Cox
          dat <- boot_backup_dat
          covariates <- tmp_covariates
          res = coxph(formula = formula(paste('Surv(OS.time, OS)~', paste(covariates,collapse=" + "))) , data = dat)
          
          #Extract results
          coef = as.data.frame(coef(summary(res)))
          ci = as.data.frame(summary(res)$conf.int)
          res = merge(coef,ci,by=0)
          res = res[,c(1:3,6,9:10)]
          
          colnames(res) = c("variable", "coefficient", "hr", "Cox.pvalue", "lower.ci", "upper.ci")
          res$hazard.ratio = paste(round(res$hr,4), " (95% CI, ", round(res$lower.ci,4), " - ", round(res$upper.ci,4), ")", sep="")
          res$Cox.pvalue = round(res$Cox.pvalue,4)
          res$prognosis = NA
          if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & res$coefficient<0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & !is.na(res$coefficient) & res$coefficient<0,]$prognosis <- "better"}
          if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & res$coefficient>0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & !is.na(res$coefficient) & res$coefficient>0,]$prognosis <- "worse"}
          res = res[,c(1,7,4,8)]
        }
      }
    } else {

      #Run multivariate Cox
      res = coxph(formula = formula(paste('Surv(OS.time, OS)~', paste(covariates,collapse=" + "))) , data = dat)
      
      #Extract results
      coef = as.data.frame(coef(summary(res)))
      ci = as.data.frame(summary(res)$conf.int)
      res = merge(coef,ci,by=0)
      res = res[,c(1:3,6,9:10)]
      
      colnames(res) = c("variable", "coefficient", "hr", "Cox.pvalue", "lower.ci", "upper.ci")
      res$hazard.ratio = paste(round(res$hr,4), " (95% CI, ", round(res$lower.ci,4), " - ", round(res$upper.ci,4), ")", sep="")
      res$Cox.pvalue = round(res$Cox.pvalue,4)
      res$prognosis = NA
      if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & res$coefficient<0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & !is.na(res$coefficient) & res$coefficient<0,]$prognosis <- "better"}
      if(nrow(res[(!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & res$coefficient>0),])>0){res[!is.na(res$Cox.pvalue) & res$Cox.pvalue<0.05 & !is.na(res$coefficient) & res$coefficient>0,]$prognosis <- "worse"}
      res = res[,c(1,7,4,8)]
    }
    
    #Merge univ_result with multi_result
    final = merge(univ_result, res, by = "variable", all = T)
    
    final[grepl("score",final$variable),]$hazard.ratio.x = as.character(logrank_result$hazard.ratio)
    final[grepl("score",final$variable),]$Cox.pvalue.x = logrank_result$log.rank.pvalue
    final[grepl("score",final$variable),]$prognosis.x = as.character(logrank_result$prognosis)

    #Replace all NAs with "----"
    final[is.na(final)] <- "----"

    #Rename columns
    colnames(final) = c("variable", "univariate.hazard.ratio", "univariate.Cox.pvalue", "univariate.prognosis",
                        "multivariate.hazard.ratio", "multivariate.Cox.pvalue", "multivariate.prognosis")

    tmp_covariates2 <- c()
    tmp_covariates3 <- c()
    covariates2 <- c("score", uni_covariates)
    new_covariates2 <- as.vector(final[[1]])

    #Match string to get from 'new covariates' the original 'covariates' name
    for(var in new_covariates2)
    {
      for(var2 in covariates2)
      {
        if(grepl(var2, var))
        {
          tmp_covariates2 <- append(tmp_covariates2, var2)
          refGroup <- substr(x = var, start = (nchar(var2) + 1), stop = nchar(var))
          tmp_covariates3 <- append(tmp_covariates3, refGroup)
        }
      }
    }

    final$variable <- tmp_covariates2
    final$reference <- tmp_covariates3
    
    final <- final[,c(1,8,2,3,4,5,6,7)]
    
    return(final)

  }, warning = function(w){ }, error = function(e){ })
}

#Create function to run multivariate Cox regression for score + clinical parameters
multiCox.model <- function(dat, univ_result, covariates){
  tmp_covariates <- c()
  #Select only relevant clinical parameters for multivariable Cox regression (p<0.2)
  covariates = c("score",covariates)

  if(roc_curve & filters)
  {
    new_covariates = c(univ_result[(univ_result$multivariate.prognosis != "----"),1])
  } else
  {
    new_covariates = c(univ_result[(univ_result$multivariate.Cox.pvalue != "----"),1])
  }
  
  #Match string to get from 'new covariates' the original 'covariates' name
  for(var in new_covariates)
  {
    for(var2 in covariates)
    {
      if(grepl(var2, var)){tmp_covariates <- append(tmp_covariates, var2)}
    }
  }
  
  tryCatch({
    
    #Run multivariate Cox
    res = coxph(formula = formula(paste('Surv(OS.time, OS)~', paste(tmp_covariates,collapse=" + "))) , data = dat)
    
    return(res)
    
  }, warning = function(w){ }, error = function(e){ })
}

#If some clinical file is provided, run multivariate analysis
if(type & clin_file != ""){

  #Print log message
  cat("Running multivariate analysis...\n\n")

  uni_covariates = colnames(clin)[4:ncol(clin)]

  #Run univariate Cox-regression for each provided clinical parameter
  cat("\tSelecting co-variables (multiple univariate analyses)...\n")
  univ_cox = suppressWarnings(univCox.test(clin,uni_covariates))
  cat("\tDone\n\n")

  #Run multivariate Cox regression if there is relevant clinical variables
  multi_cox = suppressWarnings(multiCox.test(clin, univ_cox, res_logrank, uni_covariates))

  #Write result to file
  write.table(multi_cox, paste(out,"_multiCox.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
  
  #Test proportional hazards assumptions
  cat("\tTesting proportional hazards assumptions (multivariate). Overwriting plot from 'signature score'...\n")
  multi_model = suppressWarnings(multiCox.model(dat = clin, univ_result = multi_cox, covariates = uni_covariates))
  checkPH <- test_ph_assumptions(model_object = multi_model, covariates = uni_covariates, is_multi=T)
  if (checkPH <= .05) {
    cat(paste("\tWarning: Proportional Hazards Assumptions not met (p = ", round(x = checkPH, digits = 4), "). Check plot: ",
              paste("'", out, "_ph_assumptions_plot.pdf'\n", sep=""), sep = ""))
  } else {
    cat(paste("\tProportional Hazards Assumptions met (p = ", round(x = checkPH, digits = 4), ").\n", sep = ""))
  }
  cat("\tDone\n\n")

  #Makes Forest Plot
  cat("\tMaking Forest Plot...\n")
  generate_forest_plot(model_object = multi_model, filename = paste(out, "_forest_plot.pdf", sep=""))
  cat("\tDone\n\n")
  
  #Makes bar plot if ROC curve option is TRUE
  #Check if provided files exist
  if(roc_curve2){
    cat("\tMaking BarPlot...\n")
    barplot_co_variables(plot_df = merged_table, filename = paste(out, "_frequency_bootstrap.pdf", sep = ""), covariates = uni_covariates)
    cat("\tDone\n\n")
  }
  
  #Print log message
  cat("Done\n")
  cat("\n")
}

#Print log message
cat("Analysis successfully finished\n\n")

#Calculate time elapsed to run script
end_time <- Sys.time()
elapsed_time <- difftime(time1 = end_time, time2 = start_time, units = "secs")

if (elapsed_time >= 3600) {
  cat(paste("Time to run 'survival' analysis: ", round(x = (elapsed_time[[1]] / 3600), digits = 2),
            " hours.\n", sep = ""))
} else {
  if (elapsed_time >= 60) {
    cat(paste("Time to run 'survival' analysis: ", round(x = (elapsed_time[[1]] / 60), digits = 2),
              " minutes.\n", sep = ""))
  } else {
    cat(paste("Time to run 'survival' analysis: ", round(x = elapsed_time[[1]], digits = 2),
              " seconds.\n", sep = ""))
  }
}

sink()
