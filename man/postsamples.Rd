\name{postsamples}

\alias{postsamples}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Generate samples from the posterior distribution of the GLM coefficients.
}

\description{
Generates samples from the posterior distribution of the GLM coefficients.}

\usage{
postsamples(nsample, nsubject, grid, GLMCoeffStandardized, WaveletCoefficientMatrix, 
hyperparam, pklj.bar, analysis, wave.family="DaubLeAsymm", filter.number=6, bc="periodic")
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{nsample}{Number of samples to be generated.}
  \item{nsubject}{Number of subjects included in the analysis.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is grid^2. The maximum grid value for this package is 512.}
  \item{GLMCoeffStandardized }{An array of dimension (nsubject, grid, grid), containing for each subject the standardized GLM coefficients obtained by fitting GLM to the time-series corresponding to the voxels.}
  \item{WaveletCoefficientMatrix }{A matrix of dimension (nsubject, grid^2-1), containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{hyperparam }{A vector containing the estimates of the six hyperparameters.}
  \item{pklj.bar }{A matrix of dimension (nsubject, grid^2-1), containing the piklj bar values (see References for details).}
  \item{analysis}{"MSA" or "SSA", depending on whether performing multi-subject analysis or single subject analysis.}
  \item{wave.family}{The family of wavelets to use - "DaubExPhase" or "DaubLeAsymm". Default is "DaubLeAsymm".}
  \item{filter.number}{The number of vanishing moments of the wavelet. By default 6.}
  \item{bc}{The boundary condition to use - "periodic" or "symmetric". Default is "periodic".}
}

\value{A list containing the following.
  \item{samples }{An array of dimension (nsubject,grid,grid,nsample), containing for each subject the posterior samples of the GLM coefficients.}
  \item{postdiscovery }{An array of dimension (nsubject,grid,grid), containing for each subject the posterior discovery maps of the GLM coefficients (for details see Morris et al. (2011)).}
}

\details{The wavelet computations are performed by using R package 'wavethresh'. For details, check wavethresh package help.}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.

Morris, J.S. et al. (2011). Automated analysis of quantitative image data using isomorphic functional mixed models, with
application to proteomic data. Ann. Appl. Stat. 5, 894-923.
}

\author{Nilotpal Sanyal <nisanyal@ucsd.edu>, Marco Ferreira <marf@vt.edu>
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
post.samples <- postsamples(nsample, nsubject, grid, GLMCoeffStandardized, 
WaveletCoefficientMatrix, hyperparam, pklj.bar, analysis)
dim(post.samples$samples)
#[1] 3 8 8 5
}
