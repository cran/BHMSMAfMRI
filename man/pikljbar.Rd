\name{pikljbar}

\alias{pikljbar}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute the piklj bar values of the BHMSMA model using Newton Cotes algorithm
}

\description{Computes the values of piklj bar of the BHMSMA model using Newton Cotes algorithm. For details, check References.
}

\usage{
pikljbar(nsubject, grid, WaveletCoefficientMatrix, hyperparam, analysis)
}

\arguments{
  \item{nsubject}{Number of subjects included in the analysis.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is grid^2. The maximum grid value for this package is 512.}
  \item{WaveletCoefficientMatrix }{A matrix of dimension (nsubject, grid^2-1), containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{hyperparam }{A vector containing the estimates of the six hyperparameters.}
  \item{analysis}{"MSA" or "SSA", depending on whether performing multi-subject analysis or single subject analysis.}
}

\value{A list containing the following.
  \item{pklj.bar }{A matrix of dimension (nsubject, grid^2-1), containing the piklj bar values.}
}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.
}

\author{Nilotpal Sanyal <nsanyal@stanford.edu>, Marco Ferreira <marf@vt.edu>
}

\examples{
nsubject <- 3
grid <- 8
WaveletCoefficientMatrix <- matrix(nrow=3,ncol=63)
for(i in 1:3)
 WaveletCoefficientMatrix[i,] <- rnorm(63)
hyperparam <- rep(.1,6)
analysis <- "multi"
piklj.bar <- pikljbar(nsubject, grid, WaveletCoefficientMatrix, hyperparam, analysis)
dim(piklj.bar$pklj.bar)
#[1]  3 63
}
