#Creates survival curves ##############################################

#Get command line arguments
args = commandArgs(trailingOnly=TRUE)

#Check command line argument
if (length(args) == 0) {
  stop("At least one argument must be supplied", call.=FALSE)
}

#Process command line argument
dir = args[1]

#Load packages
library(survival)
library(survminer)

#Read results from logRank tests
res = read.delim(paste(dir,"/logRank_results.txt",sep=""), header=F)

#For each line of file
for(i in 1:nrow(res)){

  #Get filename
  filename = as.character(res[i,1])

  #Get p-value
  p.val = res[i,4]

  #Get tumor type
  type = unlist(strsplit(filename,"_"))[3]

  #Define best interval
  interval = 100

  #Get gender of signature
  sig_gender = unlist(strsplit(filename,"_"))[4]

  #Get gender of samples
  if(grepl("female_data",filename)){gender="Female"}else{gender="Male"}

  #Read signature score data
  data = read.table(paste(dir,"/",filename,sep=""), header=T)

  #Get median of data signature
  sig_median = round(median(data$signature),2)

  #Transform data to categorical
  data$signature = ifelse(data$signature > median(data$signature),"high", "low")

  #If there is at least two groups to be compared, create Kaplan-meier curve
  if(nlevels(as.factor(data$signature))>1){

    pdf(paste(dir,"/plots/",gsub(".txt",".pdf",gsub("sig","plot",filename)),sep=""))
    pp = ggsurvplot(survfit(Surv(OS.time, OS) ~ signature, data=data), risk.table=TRUE, pval=paste("p =",format(p.val,scientific=T),sep=" "),font.legend=16, font.x=22, font.y=22, font.tickslab=8, pval.size=8, pval.coord=c(0,0.05), title= paste(type, sig_gender, "signature - ", gender, "samples",sep=" "), legend = c(0.7, 0.9), legend.title="", break.time.by=interval, censor=T, legend.labs = c(paste("gene signature>",sig_median, sep=""), paste("gene signature<=",sig_median, sep="")), data=data) + xlab("Survival time (days)")
    print(pp, newpage = FALSE)
    dev.off()

  }


}

