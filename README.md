
# EXDC: Estimating the Degree of Dosage Compensation at X-Linked Loci

<!-- badges: start -->
<!-- badges: end -->

The goal of EXDC is to estimate the degree of dosage compensation at X-linked loci.

## Installation

You can install the development version of EXDC like so:

``` r
if(!require(remotes)){
   install.packages("remotes")
}
remotes::install_github("yuxinyuanqt/EXDC")
```

## Example

This is a basic example which shows you how to estimate the degree of dosage compensation at X-linked loci using individual-level data:

``` r
library(EXDC)
data("Graves_Meta_analysis")
Frequen_DC(Graves_Meta_analysis$rs3827440,
           Graves_Meta_analysis$Graves_disease,
           trait_type='qualitative',
           Graves_Meta_analysis$Sex,
           method='all',con_level=0.95,
           truncation_interval=TRUE,
           truncation_lower=1,truncation_upper=2,
           description=TRUE,H0_test=TRUE,DC_0=2)
```

