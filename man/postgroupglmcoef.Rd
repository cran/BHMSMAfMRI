\name{postgroupglmcoef}

\alias{postgroupglmcoef}

\title{Obtain posterior group estimate of a 2D GLM coefficients map of a regressor
}

\description{\code{postgroupglmcoef} computes posterior group mean (or group median) of a 2D GLM coefficients map (e.g., corresponding to a single brain slice) of a regressor using the posterior means (or medians) of the corresponding wavelet coefficients from all subjects in the inverse discrete wavelet transform based on multi-subject or single subject analyses (see References).
}

\usage{
postgroupglmcoef( n, grid, glmcoefstd, postmeanwaveletcoef, 
wave.family="DaubLeAsymm", filter.number=6, bc="periodic" )
}

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
  \item{groupcoef }{A matrix of dimension (grid, grid), containing the posterior group coefficients obtained by BHMSMA methodology.}
}

\details{The wavelet transformation and reconstruction are performed by using the functions \code{\link[wavethresh]{imwd}} and \code{\link[wavethresh]{imwr}}, respectively.}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.
}

\seealso{
\code{\link{readfmridata}}, \code{\link{glmcoef}}, \code{\link{postglmcoef}}, \code{\link{postwaveletcoef}}, \code{\link{substituteWaveletCoef}}, \code{\link[wavethresh]{imwd}}, \code{\link[wavethresh]{imwr}}
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
  dim=c(n,grid^2-1))
post.groupcoef <- postgroupglmcoef(n,grid,glmcoefstd,
  postmeanwaveletcoef)
dim(post.groupcoef$groupcoef)
#[1] 8 8
}
