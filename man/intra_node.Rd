% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intra_node.R
\name{intra_node}
\alias{intra_node}
\title{Calculate the Total Weight of Self-Connections in a Network}
\usage{
intra_node(W, validate = TRUE)
}
\arguments{
\item{W}{A square, numeric matrix representing a weighted connection matrix of the network.
Self-connections are represented by the diagonal elements of the matrix.}

\item{validate}{Whether to validate the input matrix.}
}
\value{
A single numeric value representing the sum of weights of all self-connections in the network.
}
\description{
This function computes the total weight of self-connections within a network,
where a self-connection is a connection that initiates and terminates on the same node.
}
\examples{
# Create a 3x3 matrix with self-connections
W <- matrix(c(1, 2, 3,
              4, 5, 6,
              7, 8, 9), byrow = TRUE, nrow = 3)
intra_node(W)
}
