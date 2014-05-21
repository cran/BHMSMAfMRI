\name{postwaveletcoeff}

\alias{postwaveletcoeff}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Obtain posterior mean and posterior median of the wavelet coefficients using BHMSMA methodology.
}

\description{Computes posterior mean and posterior median of the wavelet coefficients using BHMSMA methodology.
}

\usage{
postwaveletcoeff(nsubject, grid, WaveletCoefficientMatrix, hyperparam, 
pklj.bar, analysis)
}

\arguments{
  \item{nsubject}{Number of subjects included in the analysis.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is grid^2. The maximum grid value for this package is 512.}
  \item{WaveletCoefficientMatrix }{A matrix of dimension (nsubject, grid^2-1), containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{hyperparam }{A vector containing the estimates of the six hyperparameters.}
  \item{pklj.bar }{A matrix of dimension (nsubject, grid^2-1), containing the piklj bar values (see Reference for details).}
  \item{analysis}{"MSA" or "SSA", depending on whether performing multi-subject analysis or single subject analysis.}
}

\value{A list containing the following.
  \item{PostMeanWaveletCoeff }{A matrix of size (nsubject, grid^2-1), containing for each subject the posterior mean of the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{PostMedianWaveletCoeff }{A matrix of size (nsubject, grid^2-1), containing for each subject the posterior median of the wavelet coefficients of all levels stacked together.}
}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.
}

\author{Nilotpal Sanyal <nsanyal_v@isical.ac.in>, Marco Ferreira <marf@vt.edu>
}

\examples{
nsubject <- 3
grid <- 8
nsample <- 5
GLMCoeffStandardized <- array(rnorm(3*8*8),dim=c(3,8,8))
WaveletCoefficientMatrix <- array(rnorm(3*63),dim=c(3,63))
hyperparam <- rep(.2,6)
pklj.bar <- array(runif(3*63),dim=c(3,63))
analysis <- "multi"
post.waveletcoeff <- postwaveletcoeff(nsubject, grid, WaveletCoefficientMatrix, 
hyperparam, pklj.bar, analysis)
dim(post.waveletcoeff$PostMeanWaveletCoeff)
#[1]  3 63
}
