\name{read.fmridata}

\alias{read.fmridata}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Read fMRI data from fMRI image files (Adopted from 'fmri').
}

\description{Reads fMRI data from fMRI image files into a 4D array.
}

\usage{
read.fmridata( directory, format, prefix, nimages, dim.image, nii=TRUE )
}

\arguments{
  \item{directory}{Location of the directory where the fMRI image files are stored. Insert within quotations ("").}
  \item{format}{The format of the data file. One of "Analyze" (.img/.hdr files), "Nifti" (.img/.hdr files or .nii files) or "Afni" (.HEAD/.BRIK files).}
  \item{prefix}{If format is "Analyze" or "Nifti", then the part of the fMRI image file name appearing before the image number. The image number is assumed to have four digit representation, that is, the lowest number is 0001 and the highest possible number is 9999. If format is "Afni", then the file name. Insert within quotations("").}
  \item{nimages}{If format is "Analyze", number of images to be read beginning from the first image. If format is "Afni", not necessary.}
  \item{dim.image}{Size of the 3D fMRI image. A vector with three elements.}
  \item{nii}{Necessary only for "Nifti" format. nii=TRUE (default) indicates the image files are in .nii files . nii=FALSE indicates the image files are .img/.hdr files.}
}

\value{A list containing the following:
  \item{fmridata }{An array of dimension (dim.image, nimages), containing the image data for all images/time-points.}
}

\details{The function uses package 'fmri' for reading from fMRI data files. For details, check fmri package help.}

\author{Nilotpal Sanyal <nisanyal@ucsd.edu>, Marco Ferreira <marf@vt.edu>
}

\seealso{
\link[AnalyzeFMRI]{f.read.analyze.volume}, \link[AnalyzeFMRI]{f.read.nifti.volume}, \link[fmri]{read.AFNI}
}

\examples{
# See package vignette.
}
