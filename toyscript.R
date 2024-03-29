#!/usr/bin/env Rscript

####Required libraries###

suppressMessages(library("TCGAbiolinks"))
suppressMessages(library("tidyverse"))

####Run in current dir####
setwd(getwd())

####Gene expression retrieval###
#Modifiable options:
query <- GDCquery(project="TCGA-GBM",                   #project-tissue option
                  data.category="Transcriptome Profiling",        #category
                  data.type="Gene Expression Quantification",     #type
                  experimental.strategy="RNA-Seq",                #experimental strategy
                  workflow.type="HTSeq - FPKM")                   #Workflow, given abundance units

GDCdownload(query, method = "api")                      #directory structure: "path"/GDCdata/TCGA-GBM/harmonized/Transcriptome_Profiling/Gen$

####Data loading and formatting###

#Gene expression#
exp_data <- GDCprepare(query)
expmatrix <- SummarizedExperiment::assay(exp_data)
expression <- t(expmatrix) %>%
   as.data.frame() %>%
   rownames_to_column(var="barcode")

#Survival status and follow_up time retrieval (clinical data)#
barcode <- exp_data$barcode
status <- exp_data$vital_status
follow_up <- exp_data$days_to_last_follow_up
tmp_clin <- cbind(barcode,status,follow_up) %>%
   as.data.frame()
tmp_clin$barcode <- as.character(tmp_clin$barcode)

###Integrating expression and clinical data###
pre_inp <- left_join(expression,tmp_clin, by="barcode") %>%
   drop_na() %>%
   filter(status!="Not Reported")

#Purge duplicated IDs from expression data
pre_inp$barcode <- substring(text = pre_inp$barcode, first=0, last=12)
pre_inp$barcode <- gsub("-","\\.", as.character(pre_inp$barcode))
dupfor<-duplicated(pre_inp$barcode)
duprev<-duplicated(pre_inp$barcode, fromLast=T)
dup <- dupfor + duprev
barfinalclin <- !(dup)
uniqindex <- which(barfinalclin)
pre_inp <- pre_inp[uniqindex,]

#Formatting and tossing genes#
len <- length(colnames(pre_inp))
aux1 <- pre_inp[,c(1,len-1,len)]

#Ensure there are no columns (genes) with only expression value equals to 0
pre_inp <- pre_inp[,2:dim(pre_inp)[2]]
exclude <- which(apply(pre_inp,2,max) < 0.001)
pre_inp <- pre_inp[,-exclude]

#Performin gene choices
len <- length(colnames(pre_inp))
set.seed(50)
random <- sample.int((len-2),50)        #Modify the number of genes as you wish (50).
aux2 <- pre_inp[,random]
input1 <- cbind(aux1,aux2)
colnames(input1) <- c("sample", "OS", "OS.time", colnames(input1[4:(dim(input1)[2])]))
input1 <- input1 %>%
   mutate(OS = if_else(OS=="Alive", 0, 1))


#############################################

###Building clinical file###

##Retrieving general info##
query <- GDCquery(project = "TCGA-GBM", data.category = "Clinical", file.type = "xml")
GDCdownload(query = query, method = "api")
clinical <- GDCprepare_clinic(query = query, clinical.info = "patient")

#Get and combine IDH and MGMT status
barcode <- exp_data$patient
barcode <- gsub("-","\\.",barcode)
IDH.status <- exp_data$paper_IDH.status
MGMT.status <- exp_data$paper_MGMT.promoter.status
first_clin <- cbind(barcode, IDH.status, MGMT.status) %>%
   as.data.frame()

#Edit and filter Drug info
clinical.drug <- GDCprepare_clinic(query = query, clinical.info = "drug")
clinical.drug$therapy_types <- ifelse(grepl(pattern = "Chemotherapy", x = clinical.drug$therapy_types, ignore.case = T), "yes", "no")
tmp_df <- clinical.drug[,c("bcr_patient_barcode", "therapy_types")]
colnames(tmp_df) <- c("barcode", "chemotherapy")
tmp_df$chemotherapy <- ifelse(tmp_df$chemotherapy == "yes", 1, 0)
final_drug <- tmp_df[order(tmp_df$barcode, -abs(tmp_df$chemotherapy) ), ] #sort by id and reverse of abs(value)
final_drug <- final_drug[!duplicated(final_drug$barcode),] # take the first row within each id
final_drug$barcode <- gsub("-","\\.", as.character(final_drug$barcode))
final_drug$chemotherapy <- ifelse(final_drug$chemotherapy == 1, "yes", "no")

#Edit and filter radiation info
clinical.radiation <- GDCprepare_clinic(query = query, clinical.info = "radiation")
clinical.radiation <- clinical.radiation[clinical.radiation$radiation_type != "",c("bcr_patient_barcode", "radiation_type")]
clinical.radiation$radiation_type <- ifelse(clinical.radiation$radiation_type != "", 1, 0)
clinical.radiation <- clinical.radiation[order(clinical.radiation$bcr_patient_barcode, -abs(clinical.radiation$radiation_type)),] #sort by i$
tmp_rad <- clinical.radiation[!duplicated(clinical.radiation$bcr_patient_barcode),] # take the first row within each id
colnames(tmp_rad) <- c("barcode", "radiotherapy")
tmp_rad$barcode <- gsub("-","\\.", as.character(tmp_rad$barcode))
tmp_rad$radiotherapy <- ifelse(tmp_rad$radiotherapy == 1, "yes", "no")

#Building final clinical toy dataset
new_clinical <- clinical[,c("bcr_patient_barcode", "history_of_neoadjuvant_treatment", "gender",
                            "age_at_initial_pathologic_diagnosis", "race_list", "ethnicity",
                            "person_neoplasm_cancer_status", "karnofsky_performance_score")]
colnames(new_clinical) <- c("barcode", "history_of_neoadjuvant_treatment", "gender",
                            "age_at_initial_pathologic_diagnosis", "race_list", "ethnicity",
                            "person_neoplasm_cancer_status", "karnofsky_performance_score")

#Ensure all variables are binary (2 categories only)
new_clinical$barcode <- gsub("-","\\.", as.character(new_clinical$barcode))
new_clinical[new_clinical == ""] <- NA
age_median <- median(new_clinical$age_at_initial_pathologic_diagnosis, na.rm = T)
new_clinical$age_at_initial_pathologic_diagnosis <- ifelse(new_clinical$age_at_initial_pathologic_diagnosis >= age_median, "HIGH", "LOW")
new_clinical$race_list <- ifelse(is.na(new_clinical$race_list), NA,
                                 ifelse(grepl(pattern = "white", x = new_clinical$race_list, ignore.case = T), "WHITE", "NO WHITE"))
karnofsky_median <- median(new_clinical$karnofsky_performance_score, na.rm = T)
new_clinical$karnofsky_performance_score <- ifelse(new_clinical$karnofsky_performance_score >= karnofsky_median, "HIGH", "LOW")

for (col in colnames(new_clinical)){
   new_clinical[[col]] <- droplevels(as.factor(new_clinical[[col]]))
}

#Merge all clinical tables
merge_therapies <- merge(tmp_rad, final_drug, by = "barcode", all = T)
outclin <- merge(new_clinical, merge_therapies, by = "barcode", all = T)
finalclin <- merge(outclin, first_clin, by = "barcode", all = T)

#Purge duplicated IDs from clinical data
dupfor <- duplicated(finalclin$barcode)
duprev <- duplicated(finalclin$barcode, fromLast=T)
dup <- dupfor + duprev
barfinalclin <- !(dup)
uniqindex <- which(barfinalclin)
finalclin <- finalclin[uniqindex,]

finalclin <- finalclin %>%
   mutate(IDH.status = if_else(IDH.status=="1", "MUT", "WT")) %>%
   mutate(MGMT.status = if_else(MGMT.status=="1", "Methylated", "Unmethylated"))

#Get patients from clinical data only if there is corresponding expression data
exp_patients <- as.factor(input1$sample)
#exp_patients <- as.vector(unique(input1$sample))
finalclin <- finalclin[finalclin$barcode %in% exp_patients,]

for (col in colnames(finalclin)){
   finalclin[[col]] <- droplevels(as.factor(finalclin[[col]]))
}

#Keep columns that have only 2 categories in clinical table
keep_cols <- c(1)
for (col in 2:ncol(finalclin))
{
   check <- nlevels(as_factor(finalclin[[col]]))
   if (check == 2)
   {
      keep_cols <- append(x = keep_cols, values = col)
   }
}
finalclin <- finalclin[,keep_cols]
colnames(finalclin) <- c("barcode", "gender", "age", "race", "ethnicity","cancer.status",
                         "karnofsky","chemotherapy","IDH.status","MGMT.status")

for (col in 2:ncol(finalclin)){
   finalclin[[col]] <- tolower(finalclin[[col]])
}

finalclin[is.na(finalclin)] <- "NA"
finalclin$chemotherapy <- ifelse(finalclin$chemotherapy == "NA", "no", ifelse(finalclin$chemotherapy == "yes", "yes", "no"))
finalclin$karnofsky <- ifelse(finalclin$karnofsky == "NA", "low", ifelse(finalclin$karnofsky == "high", "high", "low"))
finalclin[finalclin == "NA"] <- NA

##Exporting two tables: expression and clinical data##

# EXPRESSION DATA TABLE
write.table(x = input1, file = "expression.tsv", quote = F, col.names = T, row.names = F)

# CLINICAL DATA TABLE
write.table(x = finalclin, file = "clinical.tsv", col.names = T, row.names = F,
            quote = F, sep = "\t", eol = "\n", na = "NA", dec = ".")

