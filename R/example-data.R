#' Example Structural Neuroimaging and Genetic Data
#'
#' Simulated data with 632 subjects, 486 SNPs from 33 genes, 15 structural neuroimaging measures.
#'
#' @docType data
#'
#' @usage data(bgsmtr_example_data)
#'
#' @format A list with three components: "SNP_data", "SNP_groups", "BrainMeasures".
#' SNP_data is a 486-by-632 matrix containing minor allele counts for 632 subjects and 486 SNPs. 
#' SNP_groups is a vector of length 486 with labels partitioning the 486 SNPs into 33 genes.
#' BrainMeasures is a 15-by-632 matrix containing simulated volumetric and cortical thickness measures for 15 regions of interest.
#'
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(bgsmtr_example_data)
#' names(bgsmtr_example_data)
#' dim(bgsmtr_example_data$SNP_data)
#' dim(bgsmtr_example_data$BrainMeasures)
#' unique(bgsmtr_example_data$SNP_groups)
"bgsmtr_example_data"