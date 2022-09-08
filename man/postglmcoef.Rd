\name{postglmcoef}

\alias{postglmcoef}

\title{Obtain posterior estimate of a 2D GLM coefficients map of a regressor
}

\description{\code{postglmcoef} computes posterior mean (or median) of a 2D GLM coefficients map (e.g., corresponding to a single brain slice) of a regressor using the posterior mean (or median) of the corresponding wavelet coefficients in the inverse discrete wavelet transform for each subject based on multi-subject or single subject analyses (see References).
}

\usage{
postglmcoef(n, grid, glmcoefstd, postmeanwaveletcoef, 
wave.family="DaubLeAsymm", filter.number=6, bc="periodic")
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{n}{Number of subjects.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is \code{grid^2}. The maximum value of \code{grid} for this package is 512.}
  \item{glmcoefstd }{An array of dimension \code{(n,grid,grid)}, containing for each subject the standardized GLM coefficients obtained by fitting GLM to the time-series corresponding to the voxels.}
  \item{postmeanwaveletcoef }{A matrix of size \code{(n,grid^2-1)}, containing for each subject the posterior mean of the wavelet coefficients of all levels stacked together (by the increasing order of resolution level). }
  \item{wave.family}{The family of wavelets to use - "DaubExPhase" or "DaubLeAsymm". Default is "DaubLeAsymm".}
  \item{filter.number}{The number of vanishing moments of the wavelet. Default is 6.}
  \item{bc}{The boundary condition to use - "periodic" or "symmetric". Default is "periodic".}
}

\value{A list containing the following.
  \item{GLMcoefposterior }{An array of dimension \code{(n,grid,grid)}, containing for each subject the posterior means of the standardized GLM coefficients.}
}

\details{The wavelet transformation and reconstruction are performed by using the functions \code{\link[wavethresh]{imwd}} and \code{\link[wavethresh]{imwr}}, respectively.}

\seealso{
\code{\link{glmcoef}}, \code{\link{postwaveletcoef}}, \code{\link{substituteWaveletCoef}}, \code{\link[wavethresh]{imwr}}, \code{\link{postgroupglmcoef}}, \code{\link{postsamples}}
}

\author{Nilotpal Sanyal, Marco Ferreira

Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

\examples{
set.seed(1)
n <- 3
grid <- 8
glmcoefstd <- array(rnorm(n*grid*grid),
  dim=c(n,grid,grid))
postmeanwaveletcoef <- array(rnorm(n*(grid^2-1)),
  dim=c(n,(grid^2-1)))
postmeanglmcoef <- postglmcoef(n,grid,glmcoefstd,
  postmeanwaveletcoef)
dim(postmeanglmcoef$GLMcoefposterior)
#[1] 3 8 8
}
