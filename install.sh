#!usr/bin/bash

#Installation script for reboot
#Requirements:R (version >= 3.6)

Rscript -e 'install.packages(c("survival",
			       "survminer",
			       "BiocManager",
			       "optparse",
			       "OptimalCutpoints",
			       "survivalROC",
			       "forestmodel",
			       "sjstats",
			       "data.table"))'

Rscript -e 'BiocManager::install("survcomp")'

cp *.R /usr/bin
