\name{postgroupcoeff}

\alias{postgroupcoeff}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Obtain posterior group coefficients using the BHMSMA methodology.
}

\description{Computes posterior group coefficients using the BHMSMA methodology.
}

\usage{
postgroupcoeff( nsubject, grid, GLMCoeffStandardized, PostMeanWaveletCoeff, 
wave.family="DaubLeAsymm", filter.number=6, bc="periodic" )
}

\arguments{
  \item{nsubject}{Number of subjects included in the analysis.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is grid^2. The maximum grid value for this package is 512.}
  \item{GLMCoeffStandardized }{An array of dimension (nsubject, grid, grid), containing for each subject the standardized GLM coefficients obtained by fitting GLM to the time-series corresponding to the voxels.}
  \item{PostMeanWaveletCoeff }{A matrix of size (nsubject, grid^2-1), containing for each subject the posterior mean of the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{wave.family}{The family of wavelets to use - "DaubExPhase" or "DaubLeAsymm". Default is "DaubLeAsymm".}
  \item{filter.number}{The number of vanishing moments of the wavelet. By default 6.}
  \item{bc}{The boundary condition to use - "periodic" or "symmetric". Default is "periodic".}
}

\value{A list containing the following.
  \item{groupcoeff }{A matrix of dimension (grid, grid), containing the posterior group coefficients obtained by BHMSMA methodology.}
}

\details{The wavelet computations are performed by using R package 'wavethresh'. For details, check wavethresh package help.}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.
}

\author{Nilotpal Sanyal <nsanyal_v@isical.ac.in>, Marco Ferreira <marf@vt.edu>
}

\examples{
nsubject <- 3
grid <- 8
GLMCoeffStandardized <- array(rnorm(3*8*8),dim=c(3,8,8))
PostMeanWaveletCoeff <- array(rnorm(3*63),dim=c(3,63))
post.groupcoeff <- postgroupcoeff( nsubject, grid, GLMCoeffStandardized, PostMeanWaveletCoeff)
dim(post.groupcoeff$groupcoeff)
#[1] 8 8
}
