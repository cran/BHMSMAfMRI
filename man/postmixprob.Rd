\name{postmixprob}

\alias{postmixprob}

\title{Obtain estimates of the mixture probabilities defining the BHMSMA posterior wavelet coefficients distributions
}

\description{\code{postmixprob} computes the mixture probabilities (piklj.bar), which define the marginal posterior distribution of the wavelet coefficients of the BHMSMA model, using Newton Cotes algorithm for each subject based on multi-subject or single subject analyses, and returns the same (see References).
}

\usage{
postmixprob(n, grid, waveletcoefmat, hyperparam, analysis)
}

\arguments{
  \item{n}{Number of subjects.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is \code{grid^2}. The maximum value of \code{grid} for this package is 512.}
  \item{waveletcoefmat }{A matrix of dimension \code{(n,grid^2-1)}, containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{hyperparam }{A vector containing the estimates of the six hyperparameters.}
  \item{analysis}{"MSA" or "SSA", depending on whether performing multi-subject analysis or single subject analysis.}
}

\value{A list containing the following.
  \item{pkljbar }{A matrix of dimension \code{(n,grid^2-1)}, containing the piklj bar values.}
}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.
}

\author{Nilotpal Sanyal, Marco Ferreira

Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

\seealso{
\code{\link{waveletcoef}}, \code{\link{hyperparamest}}, \code{\link{postwaveletcoef}}
}

\examples{
set.seed(1)
n <- 3
grid <- 8
waveletcoefmat <- matrix(nrow=n,ncol=grid^2-1)
for(i in 1:n) waveletcoefmat[i,] <- rnorm(grid^2-1)
hyperparam <- rep(.1,6)
analysis <- "multi"
pkljbar <- postmixprob(n,grid,waveletcoefmat,hyperparam,
  analysis)
dim(pkljbar$pkljbar)
#[1]  3 63
}
