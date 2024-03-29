% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_connections.R
\name{classify_connections}
\alias{classify_connections}
\title{Classify Connections Based on Cortical and Subcortical Allocation}
\usage{
classify_connections(data_frame, cortical_allocation, subcortical_allocation)
}
\arguments{
\item{data_frame}{A data frame containing the connections to be classified.
It must include 'row' and 'col' columns representing the connections between nodes.}

\item{cortical_allocation}{A vector of nodes allocated as cortical.}

\item{subcortical_allocation}{A vector of nodes allocated as subcortical.}
}
\value{
A modified version of the input data frame with an additional column named
'connection_type' indicating the classification of each connection (Cortical-Cortical,
Subcortical-Cortical, Subcortical-Subcortical, or Unallocated).
}
\description{
This function classifies the connections in the edge data frame
produced by the function \code{process_matrices()} according to whether they are
cortical-cortical, subcortical-cortical, subcortical-subcortical, or unallocated.
It requires the data frame to contain specific columns ('row' and 'col') and
uses two vectors indicating cortical and subcortical allocations to classify the
connections.
}
\examples{
cortical_nodes <- c(1, 2, 3)
subcortical_nodes <- c(4, 5)
connections <- data.frame(row = c(1, 2, 4, 5, 3), col = c(2, 3, 5, 1, 4))
classified_connections <- classify_connections(connections, cortical_nodes, subcortical_nodes)

}
