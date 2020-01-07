##Reboot docker file##

FROM python:3.9-rc-buster

#Install R core

RUN echo "deb http://cloud.r-project.org/bin/linux/debian buster-cran35/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get -y install r-base

#Survival analysis package#

Run Rscript -e 'install.packages(c("survival","survminer","BiocManager"))'

Run Rscript -e 'BiocManager::install("survcomp")'

Run Rscript -e 'install.packages("optparse")'

Run Rscript -e 'install.packages(c("OptimalCutpoints","survivalROC","forestmodel","sjstats","data.table"))'

Run Rscript -e 'install.packages("scriptName")'

#Regression analysis package#

Run Rscript -e 'install.packages(c("penalized", "tidyverse", "hash", "R.utils"))'

#Main packages#

Run Rscript -e 'install.packages(c("argparse","extrafont"))'

#TCGA tutorial#

Run Rscript -e 'BiocManager::install("remotes")'
Run Rscript -e 'BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")'

#Tool insertion#

COPY *.R /reboot/

#Automatic running#

ENTRYPOINT ["Rscript", "/reboot/reboot.R"]
