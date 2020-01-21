#!usr/bin/bash

#Installation script for reboot
#Requirements:R (version >= 3.6)

#Packages for survival analysis
Rscript -e 'install.packages(c("survival", "survminer", "BiocManager", "optparse", "OptimalCutpoints", "survivalROC", "forestmodel", "sjlabelled", "sjstats", "sjmisc", "data.table"))'
Rscript -e 'BiocManager::install("survcomp")'

#Packages  for regression analysis 
Rscript -e 'install.packages(c("survival", "penalized", "tidyverse", "hash", "R.utils", "argparse"))'


#Packages for toy sample generation
Rscript -e 'BiocManager::install("remotes")'
Rscript -e 'BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")'

cp *.R /usr/local/bin
