##Reboot docker file##

FROM python:3.9-rc-buster

#Install R core

RUN echo "deb http://cloud.r-project.org/bin/linux/debian buster-cran35/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get -y install r-base

#Survival analysis package#

RUN Rscript -e 'install.packages(c("survival","survminer","BiocManager"))'

RUN Rscript -e 'BiocManager::install("survcomp")'

RUN Rscript -e 'install.packages("optparse")'

RUN Rscript -e 'install.packages(c("OptimalCutpoints","survivalROC","forestmodel","sjstats","data.table"))'

RUN Rscript -e 'install.packages("scriptName")'

#Regression analysis package#

RUN Rscript -e 'install.packages(c("penalized", "tidyverse", "hash", "R.utils"))'

#Main packages#

RUN Rscript -e 'install.packages(c("argparse","extrafont"))'

#TCGA tutorial#

RUN Rscript -e 'BiocManager::install("remotes")'
Run Rscript -e 'BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")'

#Tool insertion#

COPY *.R /usr/local/bin/

#Automatic running#

CMD ["Rscript", "reboot.R"]
