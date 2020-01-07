#Format input signatures (separating positive/negative coefficients) #######################################################

#Protein_coding + lncRNAs #############################################################

data = read.table("male_data2_500boot", header=T)
rownames(data) = data$g
data = data[,2,drop=F]
colnames(data) = "coefficient"
male_pos = data[data$coefficient>0,,drop=F]
male_neg = data[data$coefficient<0,,drop=F]
write.table(male_pos,"cox_signature_pc_lnc/male_pos.txt", row.names=T, col.names=T, quote=F, sep="\t")
write.table(male_neg,"cox_signature_pc_lnc/male_neg.txt", row.names=T, col.names=T, quote=F, sep="\t")

data = read.table("female_data2_500boot", header=T)
rownames(data) = data$g
data = data[,2,drop=F]
colnames(data) = "coefficient"
female_pos = data[data$coefficient>0,,drop=F]
female_neg = data[data$coefficient<0,,drop=F]
write.table(female_pos,"cox_signature_pc_lnc/female_pos.txt", row.names=T, col.names=T, quote=F, sep="\t")
write.table(female_neg,"cox_signature_pc_lnc/female_neg.txt", row.names=T, col.names=T, quote=F, sep="\t")

#Protein_coding #######################################################################

data = read.table("male_data3_500boot", header=T)
rownames(data) = data$g
data = data[,2,drop=F]
colnames(data) = "coefficient"
male_pos = data[data$coefficient>0,,drop=F]
male_neg = data[data$coefficient<0,,drop=F]
write.table(male_pos,"cox_signature_pc/male_pos.txt", row.names=T, col.names=T, quote=F, sep="\t")
write.table(male_neg,"cox_signature_pc/male_neg.txt", row.names=T, col.names=T, quote=F, sep="\t")

data = read.table("female_data3_500boot", header=T)
rownames(data) = data$g
data = data[,2,drop=F]
colnames(data) = "coefficient"
female_pos = data[data$coefficient>0,,drop=F]
female_neg = data[data$coefficient<0,,drop=F]
write.table(female_pos,"cox_signature_pc/female_pos.txt", row.names=T, col.names=T, quote=F, sep="\t")
write.table(female_neg,"cox_signature_pc/female_neg.txt", row.names=T, col.names=T, quote=F, sep="\t")
