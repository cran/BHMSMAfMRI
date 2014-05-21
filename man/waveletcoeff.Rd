\name{waveletcoeff}

\alias{waveletcoeff}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Apply discrete wavelet transform to the GLM coefficients and obtain the wavelet coefficients.
}

\description{Applies 2D discrete wavelet transform to the standardized GLM coefficient maps and returns the wavelet coefficients of all resolution levels.
}

\usage{
waveletcoeff(nsubject, grid, GLMCoeffStandardized, 
wave.family="DaubLeAsymm", filter.number=6, bc="periodic")
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{nsubject}{Number of subjects included in the analysis.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is grid^2. The maximum grid value for this package is 512.}
  \item{GLMCoeffStandardized }{An array of dimension (nsubject, grid, grid), containing for each subject the standardized GLM coefficients obtained by fitting GLM to the time-series corresponding to the voxels.}
  \item{wave.family}{The family of wavelets to use - "DaubExPhase" or "DaubLeAsymm". Default is "DaubLeAsymm".}
  \item{filter.number}{The number of vanishing moments of the wavelet. By default 6.}
  \item{bc}{The boundary condition to use - "periodic" or "symmetric". Default is "periodic".}
}

\value{A list containing the following.
  \item{WaveletCoefficientMatrix }{A matrix of dimension (nsubject, grid^2-1), containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
}

\details{The wavelet decomposition is performed by using R package 'wavethresh'. For details, check wavethresh package help.}

\author{Nilotpal Sanyal <nsanyal_v@isical.ac.in>, Marco Ferreira <marf@vt.edu>
}

\examples{
nsubject <- 3
grid <- 8
Data <- array(dim=c(3,8,8,10),rnorm(3*8*8*10))
DesignMatrix <- cbind( c(rep(c(1,0),5)), rep(1,10) )
glm.fit <- glmcoeff(nsubject, grid, Data, DesignMatrix)
GLMCoeffStandardized <- glm.fit$GLMCoeffStandardized
wavelet.coeff <- waveletcoeff(nsubject, grid, GLMCoeffStandardized)
dim(wavelet.coeff$WaveletCoefficientMatrix)
#[1]  3 63
}


