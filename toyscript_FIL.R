library("TCGAbiolinks")
library("tidyverse")

# Project can be any other TCGA tumor type such as LGG or BRCA. Same for other parameters
query <- GDCquery(project = "TCGA-GBM", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq", workflow.type = "HTSeq - FPKM")

GDCdownload(query = query, method = "api")
tran <- GDCprepare(query)

query <- GDCquery(project = "TCGA-GBM", data.category = "Clinical", file.type = "xml")
GDCdownload(query = query, method = "api")

# 'Clinical.info' can be either: drug; radiation; admin; follow_up; stage_event; new_tumor_event
clinical <- GDCprepare_clinic(query = query, clinical.info = "patient")

# Building expression file
expmatrix <- SummarizedExperiment::assay(tran)
expression <- t(expmatrix)

barcode <- clinical$bcr_patient_barcode
status <- clinical$vital_status
follow_up <- clinical$days_to_last_followup
tmp <- cbind(barcode,status,follow_up)
tmp <- tmp[complete.cases(tmp),]
tmp <- as.data.frame(tmp)

idtriple<- sapply(strsplit(rownames(expression),"-"),`[`,1:3) # 
ids=c()
for (i in seq(1,length(idtriple)-2,by=3)){
	u=i+2
	ids<-c(ids,paste(idtriple[i:u],collapse="-"))
}

rownames(expression) <- ids
tmp$barcode <- as.character(tmp$barcode)
expression <- as.data.frame(expression)
expression <- rownames_to_column(expression,var="barcode")
tmp$barcode <- gsub("-","\\.",tmp$barcode)
aux1 <- left_join(expression,tmp, by="barcode")
aux2 <- drop_na(aux1)
len <- length(colnames(aux1))
aux3 <- aux2[,c(1,len-1,len)] 
random <- sample.int(len,50)
aux4 <- aux2[,random]
aux5 <- cbind(aux3,aux4)
colnames(aux5) <- c("sample", "OS", "OS.time", colnames(aux5[4:length(colnames(aux5))]))
aux5$OS <- gsub("3","1",aux5$OS)
aux5$OS <- gsub("2","0",aux5$OS)

write.table(aux5, "expression.tsv", quote=F, col.names=T, row.names=F)

# Building clinical file
barcode <- clinical$bcr_patient_barcode
barcode <- gsub("-","\\.",barcode)
MGMT.status <-tran$paper_MGMT.promoter.status
IDH.status <- tran$paper_IDH.status

## Retrieving therapy info

# Edit and filter Drug info
clinical.drug <- GDCprepare_clinic(query = query, clinical.info = "drug")
clinical.drug$therapy_types <- ifelse(grepl(pattern = "Chemotherapy", x = clinical.drug$therapy_types, ignore.case = T), "yes", "no")
tmp_df <- clinical.drug[,c("bcr_patient_barcode", "therapy_types")]
colnames(tmp_df) <- c("barcode", "chemotherapy")
tmp_df$chemotherapy <- ifelse(tmp_df$chemotherapy == "yes", 1, 0)
final_drug <- tmp_df[order(tmp_df$barcode, -abs(tmp_df$chemotherapy) ), ] #sort by id and reverse of abs(value)
final_drug <- final_drug[!duplicated(final_drug$barcode),] # take the first row within each id
final_drug$barcode <- gsub("-","\\.", as.character(final_drug$barcode))
final_drug$chemotherapy <- ifelse(final_drug$chemotherapy == 1, "yes", "no")

# Edit and filter radiation info
clinical.radiation <- GDCprepare_clinic(query = query, clinical.info = "radiation")
clinical.radiation <- clinical.radiation[clinical.radiation$radiation_type != "",c("bcr_patient_barcode", "radiation_type")]
clinical.radiation$radiation_type <- ifelse(clinical.radiation$radiation_type != "", 1, 0)
clinical.radiation <- clinical.radiation[order(clinical.radiation$bcr_patient_barcode, -abs(clinical.radiation$radiation_type)),] #sort by id and reverse of abs(value)
tmp_rad <- clinical.radiation[!duplicated(clinical.radiation$bcr_patient_barcode),] # take the first row within each id
colnames(tmp_rad) <- c("barcode", "radiotherapy")
tmp_rad$barcode <- gsub("-","\\.", as.character(tmp_rad$barcode))
tmp_rad$radiotherapy <- ifelse(tmp_rad$radiotherapy == 1, "yes", "no")

# Building final toy dataset
new_clinical <- clinical[,c("bcr_patient_barcode", "history_of_neoadjuvant_treatment", "gender",
                            "age_at_initial_pathologic_diagnosis", "race_list", "ethnicity")]
colnames(new_clinical) <- c("barcode", "history_of_neoadjuvant_treatment", "gender",
                            "age_at_initial_pathologic_diagnosis", "race_list", "ethnicity")

new_clinical$barcode <- gsub("-","\\.", as.character(new_clinical$barcode))
new_clinical[new_clinical == ""] <- NA
age_median <- median(new_clinical$age_at_initial_pathologic_diagnosis, na.rm = T)
new_clinical$age_at_initial_pathologic_diagnosis <- ifelse(new_clinical$age_at_initial_pathologic_diagnosis >= age_median, "High", "Low")
new_clinical$race_list <- ifelse(is.na(new_clinical$race_list), NA,
                                 ifelse(grepl(pattern = "white", x = new_clinical$race_list, ignore.case = T), "white", "no white"))

for (col in colnames(new_clinical)){
  new_clinical[[col]] <- droplevels(as.factor(new_clinical[[col]]))
}

merge_therapies <- merge(tmp_rad, final_drug, by = "barcode", all = T)
outclin <- merge(new_clinical, merge_therapies, by = "barcode", all = T)

write.table(x = outclin, file = "clinical.tsv", col.names = T, row.names = F,
            quote = F, sep = "\t", eol = "\n", na = "NA", dec = ".")
