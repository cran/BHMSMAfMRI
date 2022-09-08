\name{glmcoef}

\alias{glmcoef}

\title{Fit GLM (general linear model) to the fMRI time-series of all voxels within a single 2D brain slice
}

\description{\code{glmcoef} fits a GLM to the fMRI time-series of all voxels within a single 2D brain slice for each subject, and returns standardized GLM coefficients along with their standard error for the included regressors (it does not add any intercept by itself).
}

\usage{
glmcoef(n, grid, data, designmat)
}

\arguments{
  \item{n}{Number of subjects.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is \code{grid^2}. The maximum value of \code{grid} for this package is 512.}
  \item{data}{The data in the form of an array with dimension \code{(n,grid,grid,ntime)}, where ntime is the size of the time series for each voxel.}
  \item{designmat}{The design matrix used to generate the data. An intercept column should be included unless not desired.}
}

\value{A list containing the following.
  \item{GLMCoefStandardized }{An array of dimension \code{(n,grid,grid)}, containing for each subject the standardized GLM coefficients obtained by fitting GLM to the time-series corresponding to the voxels.}
  \item{GLMCoefSE }{An array of dimension \code{(n,grid,grid)}, containing for each subject the estimated standard errors of the GLM coefficients.}
}

\references{
Friston, K.J., Holmes, A.P., Worsley, K.J., Poline, J., Frith, C.D., Frackowiak, R.S.J., 1994. Statistical parametric maps in functional imaging: a general linear approach. Hum. Brain Mapp. 2 (4), 189-210.
}

\author{Nilotpal Sanyal, Marco Ferreira

Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[stats]{lm}}, \code{\link{waveletcoef}}
}

\examples{
set.seed(1)
n <- 3
grid <- 8
ntime <- 10
designmat <- cbind(rep(1,10),c(rep(c(1,0),5)))
data <- array(dim=c(n,grid,grid,ntime),
  rnorm(n*grid*grid*ntime))
glm.fit <- glmcoef(n,grid,data,designmat)
dim(glm.fit$GLMCoefStandardized)
#[1] 3 8 8
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
%}