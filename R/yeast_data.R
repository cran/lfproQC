#' Yeast lysate - UPS1 benchmark dataset
#'
#' This dataset was given by Ramus et al., (2016). It is based on a highly complex sample (yeast lysate) spiked with different spiked amounts of the UPS1 standard mixture of 48 recombinant proteins. 
#' The original dataset contains 2644 rows of proteins and 2 groups of samples with three replicates each. 
#' The total number of missing values present in the sample is 579 (around 3.6% MVs).
#'
#' This standard proteomic dataset is suitable for benchmarking and comparing software for label-free quantification. 
#' And can also be applied to the evaluation of post-processing steps such as normalization, imputation of missing values, and statistical methods.
#'
#' Here only the portion of the dataset is taken for running the functions.
#'
#' @format A data frame with 1000 rows and 7 variables:
#' \describe{
#'   \item{Majority protein IDs}{Protein ID information}
#'   \item{A1}{1st condition, 1st technical replicate}
#'   \item{A2}{1st condition, 2nd technical replicate}
#'   \item{A3}{1st condition, 3rd technical replicate}
#'   \item{B1}{2nd condition, 1st technical replicate}
#'   \item{B2}{2nd condition, 2nd technical replicate}
#'   \item{B3}{2nd condition, 3rd technical replicate}
#'}
#' @source \doi{10.1016/j.dib.2015.11.063}
"yeast_data"
