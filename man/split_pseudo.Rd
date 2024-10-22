% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_pseudo.R
\name{split_pseudo}
\alias{split_pseudo}
\title{Split Pseudo Data}
\usage{
split_pseudo(reduced_mat, count_metadata, coi, ratio = 0.7)
}
\arguments{
\item{reduced_mat}{This takes in a traditional gene expression matrix that should have been normalized and reduced from "important_gene" function}

\item{count_metadata}{This refers to the cell type count found in metadata and should include all cell types. Should have same cell numbers as reduced_mat}

\item{coi}{This is the cell of interest, takes in the name of the cell type you are deconvoluting. Should be name from metadata}

\item{ratio}{This is the percentage split of training/test. Default to 70/30 training/testing data. Input range from 0 - 1}
}
\value{
Returns a list of training/testing data and its respective labels
}
\description{
This function prepares the training and testing data for training of the pseudo spots. The training and testing data are portions of gene
expression matrix and the labels are the portions normalized of the cell type of interest.
}
\keyword{training}
