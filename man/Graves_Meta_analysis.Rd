% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Graves_Meta_analysis.R
\docType{data}
\name{Graves_Meta_analysis}
\alias{Graves_Meta_analysis}
\title{Graves_Meta_analysis Dataset}
\format{
a data frame with 11013 rows and 4 variables:
\describe{
\item{Graves_disease}{Graves’ disease cases are coded as 1 and controls 0.}
\item{rs3827440}{for females, each genotype is coded as 0, 1 and 2 for TT, TC and CC; for males, each genotype is coded as 0 and 1 for T and C.}
\item{Sex}{the genetic sex of individuals in the sample population, where males are coded as 1 and females 2.}
\item{dataset}{there were two stages of the association analysis in Chu et al. (J Med Genet 2013, 50: 479–85), i.e., the genome-wide association study (GWAS) stage and the replication stage}
}
}
\source{
According to Chu et al. (J Med Genet 2013, 50: 479–85), SNP rs3827440 within the GPR174 gene on X chromosome was detected to be associated with the Graves’ disease.
There were two stages of the association analysis in Chu et al. (J Med Genet 2013, 50: 479–85), i.e., the genome-wide association study (GWAS) stage and the replication stage. The association between SNP rs3827440 and the Graves’ disease was identified in both of two stages and the pooled data of these two stages. There are 2941 subjects (699 males and 2242 females) in the GWAS stage and 8074 subjects (1814 males and 6260 females) in the replication stage.
According to the summary statistics in Table 1 in Chu et al. (J Med Genet 2013, 50: 479–85), we can reconstruct this dataset.
}
\usage{
Graves_Meta_analysis
}
\description{
A data frame containing X chromosome analysis for Graves’ disease
}
\references{
Chu X, Shen M, Xie F, et al. An X chromosome-wide association analysis identifies variants in GPR174 as a risk factor for Graves' disease. \emph{Journal of medical genetics}, 2013, \strong{50}: 479-485.
}
\keyword{datasets}
