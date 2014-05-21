\name{hyperparamest}

\alias{hyperparamest}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Get the estimates of the hyperparameters of the BHMSME model along with the estimate of their covariance matrix.
}

\description{Computes the MLEs of the hyperparameters of the BHMSME model following an empirical Bayes approach and the estimate of the covariance matrix of the hyperparameters.
}

\usage{
hyperparamest(nsubject, grid, WaveletCoefficientMatrix, analysis)
}

\arguments{
  \item{nsubject}{Number of subjects included in the analysis.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is grid^2. The maximum grid value for this package is 512.}
  \item{WaveletCoefficientMatrix}{A matrix of dimension (nsubject, grid^2-1), containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{analysis}{"multi" or "single", depending on whether performing multi-subject analysis or single subject analysis.}
}

\value{A list containing the following.
  \item{hyperparam }{A vector containing the estimates of the six hyperparameters of the BHMSME model.}
  \item{hyperparamVar }{Estimated covariance matrix of the hyperparameters.}
}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.
}

\author{Nilotpal Sanyal <nsanyal_v@isical.ac.in>, Marco Ferreira <marf@vt.edu>
}

\examples{
nsubject <- 3
grid <- 8
WaveletCoefficientMatrix <- array(dim=c(3,63),rnorm(3*63))
analysis <- "multi"
hyper.est <- hyperparamest(nsubject, grid, WaveletCoefficientMatrix, analysis)
}
