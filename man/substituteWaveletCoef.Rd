\name{substituteWaveletCoef}

\alias{substituteWaveletCoef}

\title{Substitute 2D wavelet transform coefficients with user-given values
}

\description{\code{substituteWaveletCoef} substitutes the wavelet coefficients stored wavelet object generated through 2D wavelet transform with user-given values and returns the modified wavelet object.
}

\usage{
substituteWaveletCoef(grid, waveletobj, values)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is \code{grid^2}. The maximum value of \code{grid} for this package is 512.}
  \item{waveletobj }{A wavelet object of class \code{imwd.object}, usually obtained by performing 2D wavelet transformion through the function \code{imwd} of package \pkg{wavethresh}.}
  \item{values}{The values with which the wavelet coefficients are to be replaced. The order should be consistent with \code{imwd.object} class.}
}

\value{A wavelet object of class \code{imwd.object} with updated wavelet coefficients.
}

\details{The maximum value of \code{grid} for this package is 512.}

\author{Nilotpal Sanyal, Marco Ferreira

Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

\seealso{
\code{\link[wavethresh]{imwd}}
}

\examples{
set.seed(1)
n <- 3
grid <- 8
ntime <- 10
designmat <- cbind( rep(1,10), c(rep(c(1,0),5)) )
data <- array(dim=c(n,grid,grid,ntime),
  rnorm(n*grid*grid*ntime))
glm.fit <- glmcoef(n, grid, data, designmat)
glmcoefstd <- glm.fit$GLMCoefStandardized[,,,1]
dwt = wavethresh::imwd(glmcoefstd[1,,],type="wavelet",
  family="DaubLeAsymm",filter.number=6,bc="periodic")
dwt

values = rnorm(grid^2-1)
dwtnew = substituteWaveletCoef(grid,dwt,values)
dwtnew
}




















