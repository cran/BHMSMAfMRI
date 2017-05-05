\name{BHMSMA}

\alias{BHMSMA}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Bayesian hierarchical multi-subject multiscale analysis of functional MRI data}

\description{Performs BHMSMA (Sanyal & Ferreira, 2012) of fMRI data using wavelet based prior that borrows strength across subjects and returns posterior smoothed versions of the fMRI data}

\usage{
BHMSMA(nsubject, grid, Data, DesignMatrix, TrueCoeff=NULL, analysis, 
wave.family="DaubLeAsymm", filter.number=6, bc="periodic")
}

\arguments{
  \item{nsubject}{Number of subjects included in the analysis.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is grid^2. The maximum grid value for this package is 512.}
  \item{Data}{The data in form of an array with dimension (nsubject,grid,grid,ntime), where ntime is the size of the time series for each voxel.}
  \item{DesignMatrix}{The design matrix used to generate the data.}
  \item{TrueCoeff}{If available, the true GLM coefficients in form of an array with dimension (nsubject,grid,grid). By default, NULL.}
  \item{analysis}{"MSA" or "SSA", depending on whether performing multi-subject analysis or single subject analysis.}
  \item{wave.family}{The family of wavelets to use - "DaubExPhase" or "DaubLeAsymm". Default is "DaubLeAsymm".}
  \item{filter.number}{The number of vanishing moments of the wavelet. By default 6.}
  \item{bc}{The boundary condition to use - "periodic" or "symmetric". Default is "periodic".}
}

\value{A list containing the following.
  \item{GLMCoeffStandardized }{An array of dimension (nsubject, grid, grid), containing for each subject the standardized GLM coefficients obtained by fitting GLM to the time-series corresponding to the voxels.}
  \item{GLMEstimatedSE }{An array of dimension (nsubject, grid, grid), containing for each subject the estimated standard errors of the GLM coefficients.}
  \item{WaveletCoefficientMatrix }{A matrix of dimension (nsubject, grid^2-1), containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{hyperparam }{A vector containing the estimates of the six hyperparameters.}
  \item{hyperparamVar }{Estimated covariance matrix of the hyperparameters.}
  \item{pklj.bar }{A matrix of dimension (nsubject, grid^2-1), containing the piklj bar values (see Reference for details).}
  \item{PostMeanWaveletCoeff }{A matrix of size (nsubject, grid^2-1), containing for each subject the posterior mean of the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{GLMcoeffposterior }{An array of dimension (nsubject, grid, grid), containing for each subject the posterior means of the standardized GLM coefficients.}
  \item{MSE }{MSE of the posterior estimates of the GLM coefficients, if the true values of the GLM coefficients are available.}
}

\details{The wavelet computations are performed by using R package 'wavethresh'. For details, check wavethresh package help.}

 \references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.
}

\author{Nilotpal Sanyal <nisanyal@ucsd.edu>, Marco Ferreira <marf@vt.edu>
}

\examples{
# Should take less than a minute to run
nsubject <- 3
grid <- 8
ntime <- 4
Data <- array(rnorm(3*8*8*4),dim=c(3,8,8,4))
DesignMatrix <- cbind(c(1,0,1,0), c(1,1,1,1))
analysis <- "multi"
BHMSMA.multi <- BHMSMA(nsubject, grid, Data, DesignMatrix, TrueCoeff=NULL, analysis)
}
