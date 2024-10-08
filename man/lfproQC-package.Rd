% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lfproQC-package.R
\docType{package}
\name{lfproQC-package}
\alias{lfproQC}
\alias{lfproQC-package}
\title{lfproQC: Quality Control for Label-Free Proteomics Expression Data}
\description{
Label-free bottom-up proteomics expression data is often affected by data heterogeneity and missing values. Normalization and missing value imputation are commonly used techniques to address these issues and make the dataset suitable for further downstream analysis. This package provides an optimal combination of normalization and imputation methods for the dataset. The package utilizes three normalization methods and three imputation methods.The statistical evaluation measures named pooled co-efficient of variance, pooled estimate of variance and pooled median absolute deviation are used for selecting the best combination of normalization and imputation method for the given dataset. The user can also visualize the results by using various plots available in this package. The user can also perform the differential expression analysis between two sample groups with the function included in this package. The chosen three normalization methods, three imputation methods and three evaluation measures were chosen for this study based on the research papers published by Välikangas et al. (2016) \doi{10.1093/bib/bbw095}, Jin et al. (2021) \doi{10.1038/s41598-021-81279-4} and Srivastava et al. (2023) \doi{10.2174/1574893618666230223150253}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/kabilansbio/lfproQC}
  \item Report bugs at \url{https://github.com/kabilansbio/lfproQC/issues}
}

}
\author{
\strong{Maintainer}: Kabilan S \email{kabilan151414@gmail.com}

Authors:
\itemize{
  \item Dr Shashi Bhushan Lal \email{sb.lal@icar.gov.in} [thesis advisor]
  \item Dr Sudhir Srivastava \email{Sudhir.Srivastava@icar.gov.in} [thesis advisor]
}

Other contributors:
\itemize{
  \item Dr Krishna Kumar Chaturvedi \email{kk.chaturvedi@icar.gov.in} [thesis advisor]
  \item Dr Yasin Jeshima K \email{yasinlab1.icar@gmail.com} [thesis advisor]
  \item Dr Ramasubramanian V \email{R.Subramanian@icar.gov.in} [thesis advisor]
  \item Dr Girish Kumar Jha \email{girish.jha@icar.gov.in} [contributor]
}

}
\keyword{internal}
