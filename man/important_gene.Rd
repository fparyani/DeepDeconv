% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/important_gene.R
\name{important_gene}
\alias{important_gene}
\title{Find Important Gene}
\usage{
important_gene(
  quant_mat,
  factor_group,
  cell_type,
  st_gene = NA,
  num_gene = 500
)
}
\arguments{
\item{quant_mat}{A gene expression matrix that has already been quantile normalized}

\item{cell_type}{This is the name of the cell whose gene signature you are looking for, should be one of the name from "groups"}

\item{st_gene}{If you are working with spatial transcriptomic data, add its gene list to ensure feasibility when running the model}

\item{num_gene}{Hyper-parameter that approximately determines number of genes to sample from all permutation of pairwise groups}

\item{groups}{The groups entered refer to the variouscell types of your dataset and are assumed to be in factored form when entered.}
}
\value{
Returns a vector of the location of the gene on the matrix inputted
}
\description{
The purpose of this function is to reduce the dimensions of a gene expression matrix by finding the most
relevant genes of a particular cell type using pairwise Wilcox test. This function assumes the distribution of your cell types in quant_mat
reflects the set of genes you are searching for.
}
\keyword{gene}
\keyword{selection}
