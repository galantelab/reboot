#!usr/local/bin/bash

#Installation script for reboot
#Requirements:R (version >= 4.0.4)

#Packages for survival analysis
Rscript -e 'install.packages(c("survival", "survminer", "BiocManager", "optparse", "OptimalCutpoints", "survivalROC", "forestmodel", "sjlabelled", "sjstats", "sjmisc", "data.table"))'
Rscript -e 'BiocManager::install("survcomp")'

#Packages  for regression analysis 
Rscript -e 'install.packages(c("penalized", "tidyverse", "hash", "R.utils", "argparse"))'
#survival package also necessary but installed previously

#Packages for toy sample generation
Rscript -e 'BiocManager::install("remotes")'
Rscript -e 'BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")'

Rscript -e 'install.packages("mice")'

cp *.R /usr/local/bin
