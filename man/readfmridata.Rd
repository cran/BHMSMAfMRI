\name{readfmridata}

\alias{readfmridata}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Import fMRI data from various fMRI image files
}

\description{\code{readfmridata} reads and imports fMRI data from various fMRI image files (Analyze, NIFTI and AFNI) into a 4D array. It is just a convenient wrapper around the data importing functions provided by the \pkg{oro.nifti} package.
}

\usage{
readfmridata( directory, format, prefix, nimages, dim.image, 
  nii=TRUE )
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

\details{The function uses package \pkg{oro.nifti} for reading from fMRI data files.}

\author{Nilotpal Sanyal, Marco Ferreira
Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

\seealso{
\code{\link[oro.nifti]{readANALYZE}}, \code{\link[oro.nifti]{readNIfTI}}, \code{\link[oro.nifti]{readAFNI}}, \code{\link{BHMSMA}}
}

\examples{
# import simmulated fMRI data from image files provided within this package
fpath <- system.file("extdata", package="BHMSMAfMRI")
untar(paste0(fpath,"/fmridata.tar"), exdir=tempdir())
data <- array(dim=c(3,32,32,9))
for(subject in 1:3)
{
  directory <- paste0(tempdir(),"/fmridata","/s0",subject,"/")
  a <- readfmridata(directory, format="Analyze", prefix=paste0("s0",subject,"_t"),
  nimages=9, dim.image=c(32,32,1))
  data[subject,,,] <- a[,,1,]
}
dim(a)
}
