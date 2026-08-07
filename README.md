# Reboot

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](...)

**Reboot** is an R package for discovering and validating prognostic molecular signatures from
high-dimensional transcriptomic data using bootstrap resampling and penalized Cox regression (LASSO).

The package is organized into two complementary modules:

* **Module I - Regression:** molecular signature discovery
* **Module II - Survival:** survival analysis for prognostic evaluation of the signature

A third function combines both modules into a complete integrated workflow.

---

## Installation

You can install the package after CRAN release:

```r
install.packages("Reboot")
```

Or, during development:

```r
devtools::install_github("galantelab/reboot")
```

---

## Overview

The package provides three main functions:

```r
rebootRegression()
rebootSurvival()
rebootComplete()
```

### Module I - `rebootRegression()`

This function:

* Reads survival + expression data
* Applies filtering steps (variance and proportional hazards)
* Performs feature selection
* Runs bootstrap LASSO Cox regression
* Returns a final molecular signature and, optionally, extra visualization resources (table, plots, and metadata)

### Module II - `rebootSurvival()`

This function:

* Applies a molecular signature to survival + expression data
* Applies filtering steps
* Calculates patient-specific scores
* Performs univariate and, optionally, multivariate survival analyses
* Returns the final model and, optionally, multiple visualization resources (tables, plots, and metadata)

### Complete workflow - `rebootComplete()`

This function sequentially executes:

1. `rebootRegression()`
2. `rebootSurvival()`

The generated molecular signature (containing `feature` and `coefficient` columns)
from module I (regression) is automatically passed to module II (survival).

---

## Input format

The expression/survival input file must contain:

* Rownames: sample IDs
* Two survival columns:

  * `OS` (event: 0/1)
  * `OS.time` (survival time)
* Remaining columns: gene or transcript expression values

Example:

```text
           OS   OS.time   Transcript1   Transcript2   Transcript3
Sample1    1    1200      5.2           3.1           0.0
Sample2    0    980       2.4           1.8           0.5
```

For multivariate analyses, an additional clinical table with matching sample IDs (rownames) is required:

```text
             Age     Stage    Gender
Sample1      Old     High     Male
Sample2      Young   Low      Female
```

Clinical variables must contain exactly two categories.

---

## Basic usage

### Module I — Regression

```r
library(Reboot)

reg_result <- rebootRegression(
  filein = example_file,
  outprefix = "my_analysis",
  bootstrap = 100,
  groupsize = 10,
  percentagefilter = 0.3,
  variancefilter = 0.01,
  followup = NULL,
  type = "transcript",
  force = TRUE,
  plots = TRUE,
  table = TRUE,
  saveJSON = TRUE,
  saveRDS = TRUE,
  report = TRUE,
  log = TRUE
)
```

### Module II — Survival

```r
surv_result <- rebootSurvival(
  filein = example_file,
  signature = reg_result,
  outprefix = "my_survival",
  multivariate = TRUE,
  clinin = clin_file,
  roc = TRUE,
  variancefilter = 0.01,
  followup = NULL,
  p.cutoff = 0.2,
  bootstrap = 100,
  force = TRUE,
  plots = TRUE,
  table = TRUE,
  saveJSON = TRUE,
  saveRDS = TRUE, 
  report = TRUE,
  log = TRUE
)
```

### Complete workflow

```r
complete_result <- rebootComplete(
  filein = example_file,
  outprefix = "my_complete_analysis",
  bootstrap = 100,
  groupsize = 10,
  percentagefilter = 0.3,
  variancefilter = 0.01,
  followup = NULL,
  type = "transcript",
  multivariate = TRUE,
  clinin = clin_file,
  roc = TRUE,
  p.cutoff = 0.2,
  force = TRUE,
  plots = TRUE,
  table = TRUE,
  saveJSON = TRUE,
  saveRDS = TRUE, 
  report = TRUE,
  log = TRUE
)
```

The parameters `bootstrap`, `variancefilter`, and `followup` are applied to both modules,
although they may operate differently during the regression and survival analysis steps.

---

## Output objects

### `rebootRegression()`

Returns an S3 object of class:

```r
reboot_signature
```

Supported methods:

```r
print()
summary()
coef()
as.data.frame()
plot()
```

---

### `rebootSurvival()`

Returns an S3 object of class:

```r
reboot_survival
```

Supported methods:

```r
print()
summary()
coef()
as.data.frame()
plot()
```

---

### `rebootComplete()`

Returns an S3 object of class:

```r
reboot_complete
```

Supported methods:

```r
print()
summary()
coef()
plot()
```

Access individual module results with:

```r
complete_result$regression
complete_result$survival
```

---

## Visualization

The package includes dedicated plotting functions for both modules. All visualization functions:

* Return editable ggplot2 or ggplot-like objects
* Support downstream customization
* Can optionally export figures to disk using `outprefix`

### Regression module

```r
reboot_histogram()
reboot_lollipop()
```

### Survival module

```r
reboot_km()
reboot_phplot()
reboot_roc()
reboot_forest()
reboot_barplot()
```

---

## Import utilities

Available import functions:

```r
read_reboot_table()
read_reboot_rds()
```

---

## Export utilities

Available export functions:

```r
write_reboot_signature()
write_reboot_score()
write_reboot_univariate()
write_reboot_multivariate()
write_reboot_metadata()
write_reboot_rds()
write_reboot_report()
```

---

## Included example datasets

The package includes toy datasets for examples and tutorials:

```r
data(toy_expression)
data(toy_clinics)
```

Additional example files are available in:

```r
system.file("extdata", package = "Reboot")
```

---

## Documentation

The package includes a complete tutorial vignette:

```r
browseVignettes("Reboot")
```

Additional resources:
* GitHub repository: https://github.com/galantelab/reboot
* Online documentation: https://galantelab.github.io/reboot/

---

## Returned S3 object structure

### `rebootRegression()`

```r
str(reg_result)
```

Contains:

* signature
* plots
* params
* metadata
* call

### `rebootSurvival()`

```r
str(surv_result)
```

Contains:

* signature
* score_cutoff
* expression_data
* clinical_data
* survival
* roc_result
* ph_result
* plots
* params
* metadata
* call

### `rebootComplete()`

```r
str(complete_result)
```

Contains:

* regression
* survival
* params
* metadata
* call

---

## Dependencies

Core package dependencies include:

* survival
* survivalROC
* survminer
* penalized
* ggplot2
* OptimalCutpoints
* sjstats

See the DESCRIPTION file for a complete and updated list of dependencies.

---

## License

GPL (>= 3)

---

## Authors

Filipe Ferreira dos Santos  
<fferreira@mochsl.org.br>

Felipe Rodolfo Camargo dos Santos  
<fefelrcs@gmail.com>

Gabriela Der Agopian Guardia  
<gguardia@mochsl.org.br>

Daniel Takatori Ohara  
<dtohara@mochsl.org.br>

Pedro Alexandre Favoretto Galante  
<pgalante@mochsl.org.br>

---

## Notes

* Runtime depends strongly on the number of bootstrap iterations and dataset size.
* Large values of `bootstrap` (e.g., 10,000) may substantially increase runtime.
* For exploratory analyses, smaller datasets and bootstrap values are recommended.
* The package is primarily designed for high-dimensional transcriptomic survival analyses.
* Users are encouraged to further investigate the REBOOT parameters and the generated signature.

---

## Future directions

Planned future extensions include:

* Interactive visualization of results
* Verbosity control
* Gain of performance with parallelization
* Explicit support to additional survival endpoints
