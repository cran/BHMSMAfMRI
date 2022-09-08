\name{postwaveletcoef}

\alias{postwaveletcoef}

\title{Obtain posterior estimates of the BHMSMA wavelet coefficients
}
\description{\code{postwaveletcoef} computes posterior mean and posterior median of the wavelet coefficients of the BHMSMA model for each subject based on multi-subject or single subject analyses (see References).
}

\usage{
postwaveletcoef(n, grid, waveletcoefmat, hyperparam, 
pkljbar, analysis)
}

\arguments{
  \item{n}{Number of subjects.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is \code{grid^2}. The maximum value of \code{grid} for this package is 512.}
  \item{waveletcoefmat }{A matrix of dimension \code{(n,grid^2-1)}, containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{hyperparam }{A vector containing the estimates of the six hyperparameters.}
  \item{pkljbar }{A matrix of dimension \code{(n,grid^2-1)}, containing the piklj bar values.}
  \item{analysis}{"MSA" or "SSA", depending on whether performing multi-subject analysis or single subject analysis.}
}

\value{A list containing the following.
  \item{PostMeanWaveletCoef }{A matrix of size \code{(n,grid^2-1)}, containing for each subject the posterior mean of the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{PostMedianWaveletCoef }{A matrix of size \code{(n,grid^2-1)}, containing for each subject the posterior median of the wavelet coefficients of all levels stacked together.}
}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.
}

\seealso{
\code{\link{waveletcoef}}, \code{\link{hyperparamest}}, \code{\link{postmixprob}}, \code{\link{postglmcoef}}
}

\author{Nilotpal Sanyal, Marco Ferreira

Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

\examples{
set.seed(1)
n <- 3
grid <- 8
nsample <- 5
waveletcoefmat <- array(rnorm(n*(grid^2-1)),
  dim=c(n,grid^2-1))
hyperparam <- rep(.2,6)
pkljbar <- array(runif(n*(grid^2-1)),
  dim=c(n,grid^2-1))
analysis <- "multi"
postwavecoef <- postwaveletcoef(n,grid,waveletcoefmat, 
hyperparam,pkljbar,analysis)
dim(postwavecoef$PostMeanWaveletCoef)
#[1]  3 63
}
