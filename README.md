
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pepquantify

<!-- badges: start -->
<!-- badges: end -->

pepquantify takes the peptide level output of the proteomics dataset and
proposes various options to pre-filter and clean up the data and prepare
for further analysis with an ‘MS-EmpiRe’ package. Also,
‘pepquantify’ can perform missing value imputation with user-defined
settings, to include peptides with a high detection rate in only one
condition, but very low in another, in the quantitative analysis.

## How to install?

You can install the development version of pepquantify from
[GitHub](https://github.com/) with:

``` r
if(!require(devtools)){
  install.packages("devtools")}

devtools::install_github("bshashikadze/pepquantify")
```

Additionally (if you do not have it already) you will need an MS-EmpiRe package installed to perform normalization & quantification
See: https://github.com/zimmerlab/MS-EmpiRe 

## Example script

This is a basic example:
1. set the working directory (folder which contains proteomics dataset)

*pepquantify read functions expect that the working directory contains proteinGroups.txt and peptides.txt (MaxQuant outputs) in case of DDA data (lfq or TMT), and main output of DIA-NN in case of the DIA analysis (to come soon). See also ?read_mqdda and ?read_diann*

### load the libraries
*The vast majority of pepquantify functions are written using tidyverse package collection, they will be imported automatically*
``` r
library(pepquantify)
library(msEmpiRe)
```


### define the function that performs data loading, normalization and quantification (MS-EmpiRe)
You need to execute this function once only in each session to have it in a global environment  

<details>
<summary>Comment</summary>

see: https://github.com/zimmerlab/MS-EmpiRe to know more about MS-EmpiRe package (also doi:10.1074/mcp.RA119.001509)  
note1: this function consists with codes which can be found in -
https://github.com/zimmerlab/MS-EmpiRe/blob/master/example.R 
note2: this is only an example code and for more information you should
refer to the documentation of an MS-EmpiRe package. 
note3: in msEmpiRe::read.standard I usually use different regext for unlisting.
This is iportant to remove unique number at the end of the protein ids which is added by pepquantify read functions and is neccessary for MS-EmpiRe read.standard function. The pattern used in the following function removes everything after the last dot (.), while the original pattern in the read.standard removes after the first dot, which is inconvinient in case of ncbi refseq protein database which has version numbers e.g. XP_123456.1
</details>


``` r
msempire_calculation <- function(data, data2 = data_raw, seed=1234, fc_threshold = 1.5) {
  
  require(magrittr)
  # read the data in the expressionset format and perform msempire normalization and quantification  
  # (https://github.com/zimmerlab/MS-EmpiRe/blob/master/example.R)
  msempiredata  <- msEmpiRe::read.standard(data[[1]], data[[2]],
                                            prot.id.generator = function(pep) unlist(strsplit(pep, "\\.[0-9]*$"))[1],
                                            signal_pattern="Intensity.*")
  
  # msempire calculations
  set.seed(seed = seed)
  msempire_results <- msempiredata  %>%
    msEmpiRe::normalize() %>%
    msEmpiRe::de.ana() %>%
    write.table(paste0(data[[3]], "/msempire_results_raw.txt"), sep = "\t", row.names = F)
  
  
  # tidy results (pepquantify package)
  pepquantify::resultstidy(data, data2,  fc_threshold = fc_threshold)}
```

### read the data
*data is loaded once*

*Keep the name as "data_raw", if you change make sure you indicate it in the next function as well (see pepquantify_funs())*

**conditions file will be generated which you should modify according to
experimental conditions. For more information please refer to the function
description (type ?read_mqdda in the R console)**

<details>
<summary>Arguments and default values</summary>

* exclude_samples:
if not empty, excludes specified sample/s from further analysis (only if necessary, e.g. after inspecting PCA)

* lfq:
if non-labelled data is loaded, lfq must be set to true, if labelling was performed (e.g. TMT) lfq should be set to false. For TMT **Reporter.intensity.corrected** is taken for quantification
</details> 

``` r
data_raw <- pepquantify::read_mqdda()
```

### prepare the dataset and perform normalization and quantification

*Two function below should be run as many times as many comparisons there are, it will generate specific folder for each comparison, e.g. if condition1 = "disease" and condition2 = "healthy", the folder will be generated automatically in the working directory named as disease_vs_healthy and all the outputs will be stored there. If you have other group(s) e.g. treated, you just copy paste two lines of code (below) and run with e.g. condition1 = "disease", condition2 = "treated" and the folder will be generated named as disease_vs_treated.*

**Also, order matters for the fold-change direction: proteins increased in abundance in the condition1 will have a positive l2fc, therefore it is prefferable that condition1 is always disease/treatment and etc. e.g. condition1 = "diabetes", condition2 = "control"**.

<details>
<summary>Arguments and default values</summary>

Options in *italics* are not (usually) necessary to change

* data:
list of two containing peptide and protein group data generated by the read functions of the pepquant package (default data_raw)

* imputation:	
if true, imputation will be performed if set to false no imputation will be performed. Generated statistics and fold-changes should be taken into account with a caution. This function is helpful to discover proteins that are missing in of the conditions while detected in another. That said it is better if imputation will be avoided in experiments with low number of samples (consider also to set second_condition to 0 (see below) in case of very small datasets) (default false)

* *n_element_peptide:*	
peptide data is the nth element of the list (change only if data is loaded manually as a list without using pepquantify read function) (default 1)

* condition1:	
name of the first condition that should be compared (note that the order matters for the fold-change direction) 

* condition2:	
name of the second condition that should be compared (note that the order matters for the fold-change direction)

* n_condition_1:	
minimum number of the valid values in the first condition (this value should be at least two, but default pepquant value is three)

* n_condition_2:	
minimum number of the valid values in the second condition (this value should be at least two, but default pepquant value is three)

* min_pep:	
minimum number of peptides for each protein (default 2)

* downshift:	
see the perseus documentation "Replace missing values from normal distribution" (default 1.8)
http://www.coxdocs.org/doku.php?id=perseus:user:activities:matrixprocessing:imputation:replacemissingfromgaussian

* width:	
see the perseus documentation "Replace missing values from normal distribution" (default 0.3)
http://www.coxdocs.org/doku.php?id=perseus:user:activities:matrixprocessing:imputation:replacemissingfromgaussian

* n_ko_like:	
minimum number of peptides that should have missing and valid value pattern (all valid in one condition, maximum 1 valid in another, or otherwise by user defined criteria (see fraction_valid and second_condition)) to be included in quantification. "ko" here does not necessarily has biological meaning, here this term is used to refer peptides that are consistently detected in one condition and not (or with extremely low rate) in another (default 2)

* fraction_valid:	
between 0-1. 1 means that imputed peptides are taken into account if they are present in all samples of one of the conditions (and max 1 in the second condition, see also option "second_condition"), 0.5 means if they are present in the half of the samples of one of the conditions. (default 1)

* second_condition:	
maximum acceptable number of valid values in other condition when fraction valid is met in the other (default 1)

* *seed:*
as values for imputation are derived randomly, seed makes sure the reproducibility (default 1234)

* fc_threshold:
minimum fold change for the protein to be considered differentially abundant (in natural scale) (default 1.5)

</details>
  
 
``` r
msempire_data <- pepquantify::pepquantify_funs(data_raw, condition1 = "name_of_condition_one", condition2 = "name_of_condition_two", imputation = FALSE)
msempire_calculation(msempire_data, fc_threshold = 1.5)
```  
  
<details>
<summary>Output</summary>

* **msempire_results_raw:**     
this is the raw results of MS-EmpiRe

* **msempire_results_tidy:**    
this is the results that has been cleaned-up and can be used for suppl tables

* **msempire_results_volcano:** 
some columns was adjusted to make it suitable for the volcano plot  

</details>  
  
