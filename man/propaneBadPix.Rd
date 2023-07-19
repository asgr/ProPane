\name{propaneBadPix}
\alias{propaneBadPix}
\alias{propanePatchPix}
\alias{propanePix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Automatically Identify Mask and Patch Bad Pixels
}
\description{
Given an image, can identify bad (cold or hot) pixels with a few parameters controlling how aggressively this is done. Either the binary mask can be returned, or the image itself with bad pixels masked, or finally the masked pixels can be replaced by good local values (patched) to fill the gaps.
}
\usage{
propaneBadPix(image, mask = NULL, inVar = NULL, smooth = 1, sigma = 10, pixcut = 1,
  cold = FALSE, hot = TRUE, dilate = FALSE, size = 3, return = "image", patch = FALSE,
  allow_write = FALSE, plot = FALSE, ...)

propanePatchPix(image, mask = NULL, smooth = 1, dilate = FALSE, size = 3,
  allow_write = FALSE, plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
An image to mask. Must be either object class Rfits_image or Rfits_pointer, i.e. as read in by \code{\link{Rfits_read_image}} or \code{\link{Rfits_point}}, or just a simple matrix.
}
  \item{mask}{
Boolean matrix or integer scalar (1,0); optional, parts of the image to mask out (i.e. ignore). This matrix mus* be the same dimensions as \option{image} where 1 means mask out and 0 means use for analysis.
}
  \item{inVar}{
Numeric matrix; the inverse-variance matrix corresponding to \option{image}. If provided this will be used to compute the sigma tension (if missing this is estimated internally). If your \option{image} has a large variance in depth then you might need to provide this (e.g. the inVar output of a \code{propaneStackWarpInVar} stack).
}
  \item{smooth}{
Numeric scalar/vector; standard deviation of the blur. This controls how sensitive we are to bad pixels and image noise. Generally should be in the range of 1-3, and for most types of bad pixel 1 is fine.
}
  \item{sigma}{
Numeric scalar/vector; the significance threshold. Higher means fewer (more extreme) pixels will be flagged. You probably want this in the range 5-20. If length 1, then the value will be used for both cold/hot, otherwise provide a vector corresponding to cold (sigma[1] / hot (sigma[2]) options.
}
  \item{pixcut}{
Numeric scalar; the maximum number of pixels that need to be grouped to be flagged as bad. Generally this should be quite low (1-5) since larger values will start to identify the cores of bright stars etc as being bad pixels. If length 1, then the value will be used for both cold/hot, otherwise provide a vector corresponding to cold (pixcut[1] / hot (pixcut[2]) options.
}
  \item{cold}{
Logical; should cold pixels be flagged? Some types of detector are prone to both cold and hot pixels (infrared) but some only experience hot pixels (hence the default is FALSE).
}
  \item{hot}{
Logical; should cold pixels be flagged? Some types of detector are prone to both cold and hot pixels (infrared) but some only experience hot pixels (hence the default is TRUE).
}
  \item{dilate}{
Logical scalar/vector; should the bad pixel mask be further dilated by a circular kernel of diameter \option{size}? If length 1, then the logic will be used for both cold/hot, otherwise provide a vector corresponding to cold (dilate[1] / hot (dilate[2]) options.
}
  \item{size}{
Numeric scalar/vector; the circular diamater of the bad pixel mask dilation (when \option{dilate} = TRUE). \option{size} must be odd, and if even then 1 will be added automatically. If length 1, then the value will be used for both cold/hot, otherwise provide a vector corresponding to cold (size[1] / hot (size[2]) options.
}
  \item{return}{
Character scalar; what to return. Options are 'image' where the input \option{image} is modified (bad pixels set to NA or patched) and returned; 'mask' where a binary mask matrix is returned, and 'loc' where the row and column locations of bad pixels are returned in a [Nbad x 2] matrix.
}
  \item{patch}{
Logical; if \option{return} = 'image', then should the bad image pixels be set to NA (FALSE) or should local average values be used to patch the information (TRUE).
}
  \item{allow_write}{
Logical; should we allow for direct writing back out to the Rfits_pointer provided to \option{image}. Only relevant if \option{image} is an Rfits_pointer and \option{return} = 'image', if both criteria are met then setting \option{allow_write} to TRUE means bad pixels will be directly modified on disk. If \option{patch} = FALSE then the bad pixels will be set to NA (so NaN on disk with CFITSIO), and if \option{patch} = TRUE then the patched values will be put there instead. Setting \option{allow_write} to TRUE is pretty aggressive, since it will actually modify data on-disk! However, if not many pixels are being masked and the source images are very large then this is much more efficient compared to writing an entire large image to disk since only the bad pixels are modified- the rest of the data on-disk is untouched.
}
  \item{plot}{
Logical; should a plot be generated of the output?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{Rwcs_image}} and \code{\link{magimage}}. Only relevant is \option{plot} = TRUE.
}
}
\details{
\code{propaneBadPix} will both identify bad pixels (cold or hot) and optionally patch them, whilst \code{propanePatchPix} will just patch in flagged (NA) pixels with local values.

In general the outputs of \code{propaneBadPix} can be used as inputs for the various stacking routines in \code{ProPane}.
}
\value{
Either a Rfits_image, Rfits_pointer or matrix depending on the combination of \option{return} and \option{allow_write}.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{propaneStackWarp}}
}
\examples{
\dontrun{
library(Rfits)
library(magicaxis)

image = Rfits_read_image(system.file('extdata/stack/image_1.fits', package="ProPane"),
  ext=2, xlo=301, xhi=500, ylo=1601, yhi=1800) #load image

set.seed(666) #set a random seed
selhot = sample(prod(dim(image)), 100)
bad = matrix(0L, dim(image)[1], dim(image)[2])
bad[selhot] = 1e3
image$imDat = image$imDat + bad

magimage(bad)
plot(image, qdiff=TRUE)

#Just make a mask:
mask_out = propaneBadPix(image, return='mask')
magimage(mask_out)

#We find all bad pixels and accidentally flag two extra ones:
magimage(bad - mask_out*1e3, qdiff=TRUE)

#And now return the original image with local patching (looks nice!):
image_fix = propaneBadPix(image, patch=TRUE, plot=TRUE, qdiff=TRUE)

#If we want our analysis to ignore some of the image we can use the mask:
mask_in = matrix(0L, 200, 200)
mask_in[1:100, 1:100] = 1L

#The bad pixels we still be present in the masked region:
propaneBadPix(image, mask_in, patch=TRUE, plot=TRUE, qdiff=TRUE)

#And just to test the direct patching:
image$imDat[selhot] = NA
plot(image, qdiff=TRUE)
propanePatchPix(image, plot=TRUE, qdiff=TRUE)

#We can further filate the bad pixels:
propanePatchPix(image, dilate=TRUE, size=5, plot=TRUE, qdiff=TRUE)
}
}

\concept{ ~mask }
