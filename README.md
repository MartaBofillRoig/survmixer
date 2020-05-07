# survmixer
Sample size and effect size calculation for overall survival based on the response rate and the survival-by-response information.


## Description

This repository contains some R functions for sample size and effect siz calculation based on a mixture distribution of responders and non-responders.

The functions included in this repository are the following:

  - *survmixture_f* to compute the survival distribution under the mixture model;
  - *survw_effectsize* to calculate the effect size according to the information on responders;
  - *survw_samplesize* to calculate the sample size according to the distributional parameters of the responders and non-responders. 


## Installation

``` r
# install.packages("devtools")
library(devtools)
install_github("MartaBofillRoig/survmixer")
```
