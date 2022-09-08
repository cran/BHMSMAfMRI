\name{postsamples}

\alias{postsamples}

\title{Obtain samples from the posterior distribution of a 2D GLM coefficient map.
}

\description{\code{postsamples} generates samples from the posterior distribution of a 2D GLM coefficient map (e.g., corresponding to a single brain slice) of a regressor in the BHMSMA model for each subject based on multi-subject or single subject analyses (see References).}

\usage{
postsamples(nsample, n, grid, glmcoefstd, waveletcoefmat, 
hyperparam, pkljbar, analysis, wave.family="DaubLeAsymm", 
filter.number=6, bc="periodic", seed)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{nsample}{Number of samples to be generated.}
  \item{n}{Number of subjects.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is \code{grid^2}. The maximum value of \code{grid} for this package is 512.}
  \item{glmcoefstd }{An array of dimension \code{(n,grid,grid)}, containing for each subject the standardized GLM coefficients obtained by fitting GLM to the time-series corresponding to the voxels.}
  \item{waveletcoefmat }{A matrix of dimension \code{(n,grid^2-1)}, containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{hyperparam }{A vector containing the estimates of the six hyperparameters.}
  \item{pkljbar }{A matrix of dimension \code{(n,grid^2-1)}, containing the piklj bar values (see References for details).}
  \item{analysis}{"MSA" or "SSA", depending on whether performing multi-subject analysis or single subject analysis.}
  \item{wave.family}{The family of wavelets to use - "DaubExPhase" or "DaubLeAsymm". Default is "DaubLeAsymm".}
  \item{filter.number}{The number of vanishing moments of the wavelet. Default is 6.}
  \item{bc}{The boundary condition to use - "periodic" or "symmetric". Default is "periodic".}
  \item{seed}{Must be a positive integer. Provide to set random number generation seed for reproducibility.}
}

\value{A list containing the following.
  \item{samples }{An array of dimension \code{(n,grid,grid,nsample)}, containing for each subject the posterior samples of the GLM coefficients.}
  \item{postdiscovery }{An array of dimension \code{(n,grid,grid)}, containing for each subject the posterior discovery maps of the GLM coefficients (for details see Morris et al. (2011)).}
}

\details{The wavelet computations are performed by using the R package \pkg{wavethresh}.}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.

Morris, J.S. et al. (2011). Automated analysis of quantitative image data using isomorphic functional mixed models, with application to proteomic data. Ann. Appl. Stat. 5, 894-923.
}

\author{Nilotpal Sanyal, Marco Ferreira

Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

\seealso{
\code{\link{readfmridata}}, \code{\link{glmcoef}}, \code{\link{waveletcoef}}, \code{\link{hyperparamest}}, \code{\link{postmixprob}}, \code{\link{postwaveletcoef}}, \code{\link{substituteWaveletCoef}}, \code{\link{postglmcoef}}, \code{\link[wavethresh]{imwd}}, \code{\link[wavethresh]{imwr}}
} 

\examples{
set.seed(1)
n <- 3
grid <- 8
nsample <- 5
glmcoefstd <- array(rnorm(n*grid*grid),
  dim=c(n,grid,grid))
waveletcoefmat <- array(rnorm(n*(grid^2-1)),
  dim=c(n,(grid^2-1)))
hyperparam <- rep(.2,6)
pkljbar <- array(runif(n*(grid^2-1)),
  dim=c(n,(grid^2-1)))
analysis <- "multi"
postsample <- postsamples(nsample,n,grid,glmcoefstd, 
waveletcoefmat, hyperparam,pkljbar,analysis,seed=1)
dim(postsample$samples)
#[1] 3 8 8 5
}
