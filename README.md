# survmixer

Sample size and effect size calculation for overall survival given the information on short-term binary endpoints.

In this work, we consider a mixture model that relates both the survival
and binary endpoints.
We distinguish between patients who respond to the binary endpoint, called responders, ant those who don't, non-responders. 
We use the difference of the restricted mean survival times (RMSTs) as the basis of the comparison between arms. 

The expected effect size (RMSTs difference) and sample size are then calculated on the basis of the response rate of the binary endpoint as well as on the survival functions for responders and non-responders in each treatment arm.


## Description

This repository contains R functions for sample size and effect size calculation according to different set of parameters.

The functions included in this repository are the following:

  - *survmixture_f* to compute the survival distribution under the mixture model;
  - *survm_effectsize* to calculate the effect size (in terms of the RMST difference) according to the information on responders;
  - *survm_samplesize* to calculate the sample size according to the distributional parameters of the responders and non-responders. 


## Installation

``` r
# install.packages("devtools")
library(devtools)
install_github("MartaBofillRoig/survmixer")
```

## References

This repository also contains the source files of the preprint:

- "Design of phase III trials with long-term survival outcomes based on short-term binary results". Marta Bofill Roig, Yu Shen, Guadalupe Gómez Melis. (2020). 
https://arxiv.org/abs/2008.12887


