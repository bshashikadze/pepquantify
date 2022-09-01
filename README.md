
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pepquantify

<!-- badges: start -->
<!-- badges: end -->

pepquantify Takes the peptide level output of the proteomics dataset and
proposes various options to pre-filter and clean up the data and prepare
for the further analysis with an ‘MS-EmpiRe’ package. Also,
‘pepquantify’ can perform missing value imputation with user-defined
settings, to include peptides with the high detection rate in only one
condition, but low in another, in the quantitative analysis.

## Installation

You can install the development version of pepquantify from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bshashikadze/pepquantify")
```

## Example

This is a basic example:

### load the libraries

``` r
library(pepquantify)
library(msEmpiRe)
```

### define the function that performs data loading, normalization and quantification (MS-EmpiRe)

see: <https://github.com/zimmerlab/MS-EmpiRe> note1: this function
consists with codes which can be found in -
<https://github.com/zimmerlab/MS-EmpiRe/blob/master/example.R> note2:
that this is only an example code and for more information you should
refer to the documentation of an MS-EmpiRe package.

``` r
msempire_calculation <- function(data, data2 = data_raw, seed=1234, fc_threshold = 1.5) {
  
  require(magrittr)
  
  # read the data in the expressionset format and perform msempire normalization and quantification  
    # (https://github.com/zimmerlab/MS-EmpiRe/blob/master/example.R)
  msempire_data  <- msEmpiRe::read.standard(msempire_data[[1]],msempire_data[[2]],
                                            prot.id.generator = function(pep) unlist(strsplit(pep, "\\.[0-9]*$"))[1],
                                            signal_pattern="Intensity.*")
  
  # msempire calculations
  set.seed(seed = seed)
  msempire_results <- msempire_data  %>%
    msEmpiRe::normalize() %>%
    msEmpiRe::de.ana() %>%
    write.table(paste0(data[[3]], "/msempire_results_raw.txt"), sep = "\t", row.names = F)
  
  # tidy results
  pepquantify::resultstidy(data, data2,  fc_threshold = fc_threshold)
}
```

### read the data

conditions file will be generated which you should modify according to
experimental conditions for more information please read the function
description by ?read_mqdda

Arguments:
exclude_samples
if not empty, excludes specified sample/s from further analysis (only if necessary, e.g. after inspecting PCA)

lfq	
if non-labelled data is loaded, lfq must be set to true if labelling was performed (e.g. TMT) lfq should be set to false. For TMT Reporter.intensity.corrected is taken for quantification


``` r
data_raw <- pepquantify::read_mqdda()
```

### prepare the dataset and perform quantification
'imputation' true/false
``` r
msempire_data <- pepquantify_funs(data_raw, condition1 = "name_of_condition_one", condition2 = "name_of_condition_two")
msempire_calculation(msempire_data, fc_threshold = 1.5)
```
