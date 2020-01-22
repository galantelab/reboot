## Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Dependencies](#dependencies)
4. [Usage and options](#usage-and-options)
5. [Inputs](#inputs)
6. [Toy example](#toy-example)
7. [Outputs](#outputs)

## Overview

Reboot is a modular tool, developed in R 3.6, comprising two main functionalities: regression and survival.
It was built to provide the freedom of choice for regression and survival analysis; only regression analysis; or only survival analysis	

![](Paper_figure.png)
*Figure 1: Reboot worflow. First module runs a regression analysis. Second module runs survival analysis based on a generated signature* 

Regression module produces significance coefficients (genetic signature) based on a Cox regression for genes and/or transcripts. Dimension is a critical variable for algorithm convergence in regression procedures. 
{: .alert .alert-info .text-justify}

The survival module produces and applies a score for all individuals based on the performed regression module. Both univariate and multivariate survival analyses use by default the median score value as cutoff for stratification of patients in high and low score genetic signature, unless the ROC option is chosen. In that case, the cutoff value is based on the ROC curve using NNE (Nearest Neighbour Estimate) method and the Youden statistics, where J = [sensitivity + (specificity -1)]. If more than one J coefficients is present, then the first one is chosen.
{: .alert .alert-info .text-justify}

Reboot also offers the multivariate option, where other clinical variables such as therapy, age, gender, among others can be included for a multivariate survival analysis. Multiple univariate analysis are made and only variables with a p-value <= 0.2 are selected for the final multivariate analysis. Statistical tests are performed in order to evaluate the relevance of each co-variable as a prognostic factor of a given event (overall / progression-free / recurrence-free survival).
{: .alert .alert-info .text-justify}

Additionally, if the ROC option is chosen along with the multivariate option, the multivariate analysis is done with a bootstrap resampling method if the clinical dataset provided passes the filters: (i) final dataset with at least 70% of the original one (NAs filter) and; (ii) the frequency of the less abundant category for each co-variable is not less than 20% (proportion filter). Otherwise, a multivariate analysis is performed without the bootstrap method. After 100 iterations, the relevance frequency of each co-variable with the event is calculated. Several plots are drawn for variables whose frequencies are at least 25% 
{: .alert .alert-info .text-justify}

## Installation

Reboot can be obtained from github, through direct installation or docker container (recommended).

1. Via docker
	
	This method works on any distribution or operational systems, as long as [docker](https://docs.docker.com/install){:target="_blank"} is installed.

	```docker pull galantelab/reboot```

2. Direct

	This method requires previous installation of  [R version >= 3.6](https://cran.r-project.org/doc/manuals/r-devel/R-admin.html){:target="_blank"}:

	```git clone https://github.com/galantelab/reboot.git```


	```sudo sh reboot/install.sh```


## Dependencies

In order to work properly, the following packages are necessary (included in the installation procedures):

* survival
* survminer
* [BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html){:target="_blank"}
* survcomp
* optparse
* OptimalCutpoints
* survivalROC
* forestmodel
* sjlabelled
* sjmisc
* sjstats
* data.table
* penalized
* tidyverse
* hash
* R.utils
* argparse
* extrafont
* remotes (from BiocManager)
* BioinformaticsFMRP/TCGAbiolinks (from BiocManager)

## Commands and options 

Reboot works with a command/subcommands structure:

```reboot.R [subcommand] <options>```


Subcommands may be invoked by the help menu:

```docker run --rm galantelab/reboot reboot.R -h```  or optionally:

```reboot.R -h```


In summary, 3 subcommands are available:


| **regression** | generates signature through multivariate cox regression analysis |
| **survival** | applies signature in survival analysis |
| **complete** | generates and applies signature | 

<br>

1. **regression**

   In order to generate a genetic signature, run the following: 


   ```docker run --rm galantelab/reboot regression <options>``` , or optionally:


   ```reboot.R regression <options>```


   Regression options are: 

   | short option | full option | Description|
   | ----------------------- | ----------------------------- | ----------------------- |
   | -I | -\-filein | Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival paramenters|
   | -O | [--]outprefix |  Output file prefix. Default: reboot |
   | -B | [--]bootstrap | Number of iterations for bootstrap simulation (int). Default: 1 |
   | -G | [--]groupsize | Number of genes/transcripts to be selected in each bootstrap simulation (int). Default: 3 |
   | -P | [--]pcentfilter | Percentage of correlated gene/transcript pairs allowed in each iteration. Default: 0.3 |
   | -V | [--]varfilter | Minimum normalized variance (0-1) required for each gene/transcript among samples (double). Default: 0.01 |
   | -h | [--]help      | Show this help message and exit |

   <br>	

2. Survival option

   In order to generate a survival analysis, run the follwing:


   ```docker run --rm galantelab/reboot reboot.R survival <options>```, or optionally:


   ```reboot.R survival <options>```


   If you are using docker container:


   Survival options are:

   | short option | full option | Description|
   | ------------------------ | ----------------------------- | ----------------------- |
   | -I | - -filein | Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival paramenters|
   | -O | - -outprefix |  Output file prefix. Default: reboot |
   | -M | - -multivariate | If clinical variables should be included, choose -M. This option is tied with -C option |
   | -C | - -clinical | Tab separated values (tsv) file containing binary categorical variables only. Required if -M option is chosen |
   | -R | - -roc | If continuous variables should be categorized according to a ROC curve instead of median, choose -R |
   | -S | - -signature | Tab separated values (tsv) file containing a set of genes/transcripts and corresponding cox coefficients |
   | -h | [--]help      | Show this help message and exit |

   <br>

4. Complete option

   All sub-options are provided by: (give a more detailed table after command)

   ```reboot.R complete -h```

   If you are using docker container:

   ```docker run --rm  galantelab/reboot reboot.R complete -h```

   In summary, options are:

   | short option | full option | Description|
   | ----------------------- | ----------------------------- | ----------------------- |
   | -I | - -filein | Input file name. Tab separated values (tsv) file containing genes/transcripts expression and survival paramenters|
   | -O | - -outprefix |  Output file prefix. Default: reboot |
   | -B | - -bootstrap | Number of iterations for bootstrap simulation (int). Default: 1 |
   | -G | - -groupsize | Number of genes/transcripts to be selected in each bootstrap simulation (int). Default: 3 |
   | -P | - -pcentfilter | Percentage of correlated gene/transcript pairs allowed in each iteration. Default: 0.3 |
   | -V | - -varfilter | Minimum normalized variance (0-1) required for each gene/transcript among samples (double). Default: 0.01 |
   | -M | - -multivariate | If clinical variables should be included, choose -M. This option is tied with -C option |
   | -C | - -clinical | Tab separated values (tsv) file containing binary categorical variables only. Required if -M option is chosen |
   | -R | - -roc | If continuous variables should be categorized according to a ROC curve instead of median, choose -R |

   <br>

## Inputs
		
1. It is required quantitative information of atributes, including status and follow up time which must be in tsv format. This file is required for both regression and survival analysis:

   | Sample ID | OS | OS.time | Feature1 | Feature2 | ... |   
   |---|---|---|---|---|---|
   | Sample1 | OS status | OS.time value | Feature1 value | Feature2 value | ... | 
   | Sample2 | OS status | OS.time value | Feature1 value | Feature2 value | ... |
   | ... | ... | ... | ... | ... | ... |

   <br>

2. In case survival option is invoked, clinical information is also necessary in a tsv format, in the following configuration (all clinical variables MUST be categorical and present ONLY 2 classes):

   | Sample ID | clinical_variable1 | clinical_variable2 | clinical_variable3 | ... |   
   |---|---|---|---|---|
   | Sample1 | clin_var1(status) | clin_var2(status) | clin_var3(status) | ... | 
   | Sample2 | clin_var1(status) | clin_var2(status) | clin_var3(status) | ... |
   | ... | ... | ... | ... | ... |

   <br>

## Toy example

   In order to ilustrate usage, a toy script is provided to download and format expression and clinical data of glioblastoma patients from TCGA.
   Running the following code in the reboot directory provides both inputs:
	
   ```toyscript.R```

   It is also possible to obtain toy datasets from reboot docker image, using the following:

   ```docker run --env MYID=$(id -u) --rm -ti -v $(pwd):$(pwd) -w $(pwd) galantelab/reboot toyfordocker.R```

   This command returns 2 tsv files, mentioned above, called expression.tsv and clinical.tsv. A MANIFEST.txt file and a set of expression and clinical data are also created, as intermediates of TCGA dowload process.
   The composition of expression dataset comprises clinical variables: OS (survival status) and OS.time (follow up time) and 50 random picked gene expression (FPKM).

   Finally, reboot can be run in the complete mode:

   ```Rscript reboot.R complete -I expression.tsv -O toy -B 100 -G 10 -M -C clinical.tsv -R```

   Docker:
  
   ```docker run -u $(id -u):$(id -g) --rm -v $(pwd):$(pwd) -w $(pwd) galantelab/reboot reboot.R complete -I expression.tsv -O toy -B 100 -G 10 -M -C clinical.tsv -R``` 
    
## Outputs
 
1. Regression option

   Generated outputs comprise one log file, a tsv containing relevance coefficients and 2 plots. The text file is in the following format:
	
   | Feature name | coefficient | 
   | --- | --- |
   | Feature1 | coefficient1 | 
   | Feature2 | coefficient2 |
   | ... | ... |

   The plots in the output are a histogram with the distribution of the coefficients and a lollipop plot with the most relevant coefficients.


2. Survival option

   Survival analysis may be done in multivariate or univariate mode.

   #### Univariate mode

      If the analysis is performed in univariate mode, reboot returns a log and a lograng.txt file, containing the statistics of the signature:

      | feature | coefficient | hazard.ratio | log.rank.pvalue | low.high.samples | median.survival.low | median.survival.low | prognosis |
      | --------------- | ----------- | ------------ | --------------- | ---------------- | ------------------- | ------------------- | --------- | 
      | signature_score | coefficient value | hazard.ratio value | log.rank.pvalue value | low.high.samples value | median.survival.low value | median.survival.low value | prognosis value |
       
      Plots returned for this mode are proportional hazard assumptions plot, containing a Schoenfeld test and a Kaplan Mayer plot.

   #### Multivariate mode
 
      If the analysis is performed in multivariate mode, reboot returns a multicox.txt file containing the statistics of the signature and all other clinical variables, besides the log file:

      | feature | coefficient | hazard.ratio | log.rank.pvalue | low.high.samples | median.survival.low | median.survival.low | prognosis |
      | --------------- | ----------- | ------------ | --------------- | ---------------- | ------------------- | ------------------- | --------- |
      | signature_score | coefficient value | hazard.ratio value | log.rank.pvalue value | low.high.samples value | median.survival.low value | median.survival.low value | prognosis value |
      | clin variable 1 | coefficient value | hazard.ratio value | log.rank.pvalue value | low.high.samples value | median.survival.low value | median.survival.low value | prognosis value |
      | clin variable 2 | coefficient value | hazard.ratio value | log.rank.pvalue value | low.high.samples value | median.survival.low value | median.survival.low value | prognosis value|
      | ... | ... | ... | ... | ... | ... | ... | ... |

      Plots returned for this mode are: a forest plot for all clinical variables considered, a Kaplan Meier plot and a proportional hazard assumption plot metioned in the univariate mode. If the option --ROC is also selected then the split of the patients in high and low signature score is done based on a ROC curve, also provided. 

3. Complete option

   This option is the integration of regression analysis with survival analysis. Thus, all plots and tables generated by regression and survival analysis will be produced according to the options desired. A single log file is created for this option
