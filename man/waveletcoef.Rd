\name{waveletcoef}

\alias{waveletcoef}

\title{Apply discrete wavelet transform (DWT) to a 2D GLM coefficient map of a regressor
}

\description{\code{waveletcoef} applies DWT to a 2D GLM coefficient map (e.g., corresponding to a single brain slice) of a regressor for each subject, and returns the wavelet coefficients at all resolution levels. This function wraps around the wavelet transformation function \code{\link[wavethresh]{imwd}} of the \pkg{wavethresh} package. 
}

\usage{
waveletcoef(n, grid, glmcoefstd, wave.family="DaubLeAsymm", 
  filter.number=6, bc="periodic")
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{n}{Number of subjects.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is \code{grid^2}. The maximum value of \code{grid} for this package is 512.}
  \item{glmcoefstd }{An array of dimension \code{(n,grid,grid)}, containing for each subject the standardized GLM coefficients of a regressor obtained by fitting GLM to the time-series corresponding to the voxels.}
  \item{wave.family}{The family of wavelets to use - "DaubExPhase" or "DaubLeAsymm". Default is "DaubLeAsymm".}
  \item{filter.number}{The number of vanishing moments of the wavelet. Default is 6.}
  \item{bc}{The boundary condition to use - "periodic" or "symmetric". Default is "periodic".}
}

\value{A list containing the following.
  \item{WaveletCoefficientMatrix }{A matrix of dimension \code{(n,grid^2-1)}, containing for each subject the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
}

\details{The wavelet decomposition is performed by using the function \code{\link[wavethresh]{imwd}}.}

\author{Nilotpal Sanyal, Marco Ferreira

Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

\seealso{
\code{\link[wavethresh]{imwd}}, \code{\link{hyperparamest}}
}

\examples{
set.seed(1)
n <- 3
grid <- 8
ntime <- 10
designmat <- cbind( rep(1,10), c(rep(c(1,0),5)) )
data <- array(dim=c(n,grid,grid,ntime),
  rnorm(n*grid*grid*ntime))
glm.fit <- glmcoef(n,grid,data,designmat)
glmcoefstd <- glm.fit$GLMCoefStandardized[,,,1]
wavecoef <- waveletcoef(n,grid,glmcoefstd)
dim(wavecoef$WaveletCoefficientMatrix)
#[1]  3 63
}


