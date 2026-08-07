##### PREPARE R ENV #####
# Installs REBOOT R package
# install.packages("../Reboot_1.2.0.tar.gz", repos = NULL, type = "source")

# Loads REBOOT R package
library(Reboot)
#########################

##### LOADS TOY DATASETS #####
# Loads internal dataset
# Alternatively, users can run read_reboot_table() --> see vignette (end of code) for details
data("toy_expression")
data("toy_clinics")

# Users have the option to load previous results saved in RDS
# result_complete <- read_reboot_rds(file = "extra_toy_testComplete.rds")
##############################

##### EXECUTES REBOOT #####
# Runs REBOOT complete workflow --> modules I (regression) + II (survival)
result_complete <- rebootComplete(
  filein = toy_expression,
  outprefix = "toy_testComplete",
  bootstrap = 10,
  groupsize = 10,
  percentagefilter = 0.3,
  variancefilter = 0.01,
  followup = NULL,
  type = "transcript",
  multivariate = T,
  clinin = toy_clinics,
  roc = T,
  p.cutoff = 0.2,
  force = T,
  plots = T,
  table = T,
  saveJSON = T,
  saveRDS = T,
  report = T,
  log = T
)
###########################

##### INSPECT OBJECTS #####
# Explores the "reboot_complete" S3 object
result_complete$regression  # "reboot_signature" S3 object
result_complete$survival  # "reboot_survival" S3 object
result_complete$params
result_complete$metadata
result_complete$call

# Tests S3 methods
print(result_complete)
print(summary(result_complete))
summary(result_complete)
coef(result_complete)
plot(result_complete)
###########################

##### FIGURES & TABLES #####

### MODULE I - REGRESSION

# Exports plots manually
reboot_histogram(result_complete$regression$signature , "extra_toy_testComplete")
reboot_lollipop(result_complete$regression$signature, "extra_toy_testComplete")

# Exports signature manually
write_reboot_signature(result_complete$regression, "extra_toy_testComplete_signature.txt")

### MODULE II - SURVIVAL

# Exports main plots manually
reboot_roc(x = result_complete$survival$roc_result, outprefix = "extra_toy_testComplete")
reboot_km(data = result_complete$survival$expression_data, cutoff = result_complete$survival$score_cutoff,
          pval = result_complete$survival$survival$univariate$table$log.rank.pvalue, outprefix = "extra_toy_testComplete")
reboot_phplot(x = result_complete$survival$ph_result, is_multi = TRUE, outprefix = "extra_toy_testComplete")
reboot_forest(model_object = result_complete$survival$survival$multivariate$model, outprefix = "extra_toy_testComplete")

# Exports updated signature manually
write_reboot_signature(result_complete$survival, "extra_toy_testComplete_signature_updated.tsv")

# Exports scores (with or without clinical data) manually
write_reboot_score(result_complete$survival, "extra_toy_testComplete_scoreCont.tsv")
write_reboot_score(result_complete$survival, "extra_toy_testComplete_scoreCat.tsv", clinics = TRUE)

# Exports uni and multivariate survival results manually
write_reboot_univariate(result_complete$survival, "extra_toy_testComplete_uniCox.txt")
write_reboot_multivariate(result_complete$survival, "extra_toy_testComplete_multiCox.txt")

### EXTRA

# Exports metadata manually for reproducibility
write_reboot_metadata(module = "rebootComplete", outprefix = "extra_toy_testComplete", command = "none",
                      parameters = list(filein = "data.frame", clinin = "data.frame", signature = "data.frame",
                                        outprefix = "extra_toy_testComplete", multivariate = T, roc = T, groupsize = 10,
                                        percentagefilter = 0.3, variancefilter = 0.01, followup = NULL, p.cutoff = 0.2,
                                        bootstrap = 10, type = "transcript", force = T, plots = T,
                                        table = T, saveJSON = T, saveRDS = T, report = T, log = T))

# Saves result in RDS for flexibility
write_reboot_rds(object = result_complete, file = "extra_toy_testComplete.rds")

# Generates final report for collaborators
write_reboot_report(object = result_complete, format = "PDF", file = "extra_toy_testComplete_report_complete.pdf")
############################

##### VIGNETTE #####
# Explores the REBOOT tutorial
browseVignettes("Reboot")
vignette("reboot-tutorial", package = "Reboot")
####################
