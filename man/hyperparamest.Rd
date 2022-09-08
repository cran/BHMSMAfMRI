\name{hyperparamest}

\alias{hyperparamest}

\title{Obtain estimates of the hyperparameters of the BHMSMA model
}

\description{\code{hyperparamest} computes the MLEs (maximum likelihood estimates) of the hyperparameters of the BHMSMA model using an empirical Bayes approach for multi-subject or single subject analyses, and returns the hyperparameters estimates along with their covariance matrix estimate (see References).
}

\usage{
hyperparamest(n, grid, waveletcoefmat, analysis)
}

\arguments{
  \item{n}{Number of subjects.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is \code{grid^2}. The maximum value of \code{grid} for this package is 512.}
  \item{waveletcoefmat}{A matrix of dimension \code{(n,grid^2-1)}, containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{analysis}{"multi" or "single", depending on whether performing multi-subject analysis or single subject analysis.}
}

\value{A list containing the following.
  \item{hyperparam }{A vector containing the estimates of the six hyperparameters of the BHMSMA model.}
  \item{hyperparamVar }{Estimated covariance matrix of the hyperparameters.}
}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.
}

\author{Nilotpal Sanyal, Marco Ferreira

Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

\seealso{
\code{\link{waveletcoef}}, \code{\link[stats]{nlminb}}, \code{\link{postmixprob}}
}

\examples{
set.seed(1)
n <- 3
grid <- 8
waveletcoefmat <- array(dim=c(n,grid^2-1),
  rnorm(n*(grid^2-1)))
analysis <- "multi"
hyperest <- hyperparamest(n,grid,waveletcoefmat,analysis)
hyperest$hyperparam
# [1]  1.00000  1.00000  1.00000  1.00000  0.00000 28.37678
}
