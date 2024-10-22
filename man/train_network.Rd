% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_network.R
\name{train_network}
\alias{train_network}
\title{Train Neural Network}
\usage{
train_network(train_cell_data, train_labels)
}
\arguments{
\item{train_cell_data}{Transpose of regular gene expression matrix of data that requires training. split_pseudo does the transposing automatically}

\item{train_labels}{This is the training data which is the portion of cell type to predict.}
}
\value{
Returns a trained network. Can be saved as an .h5 file using TensorFlow
}
\description{
The purpose of this function is to traing the pseudo spot data on a neural network. This function uses an external depedency for setting
up the neural network: Keras. Can freely change parameters of neural network. Feel free to refer to Keras guide on NN.
}
\keyword{training}
