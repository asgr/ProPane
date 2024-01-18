\name{propaneStackFlat}
\alias{propaneStackFlat}
\alias{propaneStackFlatInVar}
\alias{propaneStackFlatFunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Flat Image Stacking
}
\description{
Stacks multiple aligned (flat) images.
}
\usage{
propaneStackFlatInVar(image_list = NULL, sky_list = NULL, skyRMS_list = NULL,
  mask_list = NULL, magzero_in = 0, magzero_out = 23.9, masking = 'and')

propaneStackFlatFunc(image_list = NULL, imager_func = NULL, na.rm = TRUE, weights = NULL,
  increasing = TRUE, probs = c(0.159, 0.5, 0.841), ondisk = FALSE,
  cores = floor(detectCores()/2), multitype = 'fork', chunk = 1e3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_list}{
List; required. Each list element is a numeric matrix, Rfits_image or Rfits_pointer representing the images to be stacked. Note, image NAs are treated as masked pixels.
}
  \item{sky_list}{
List; optional. Each list element is a numeric matrix representing the sky to be subtracted.
}
  \item{skyRMS_list}{
List; optional. Each list element is a numeric matrix representing the sky-RMS to weight the stack with.
}
  \item{mask_list}{
List optional. Each list element is a integer matrix representing a mask for each image in \option{image_list}, where 1 means mask out and 0 means use for analysis. This can be used in conjuction with setting pixels in \option{image_list} to NA.
}
  \item{magzero_in}{
Numeric vector; the input mag-zero points. If length 1 then it is assumed all input frames have the same mag-zero point.
}
  \item{magzero_out}{
Numeric scalar; the output mag-zero point desired.
}
  \item{masking}{
Character scalar; what to do with masked pixels (NAs in the \option{image}). If 'or'/'|'/'||' a pixel is masked if *any* of the images being stacked have a masked pixel (NA in the \option{image_list}) at that location, if 'and'/'&'/'&&' then a pixel is masked if *all* of the images being stacked have a masked pixel at that location.
}
  \item{imager_func}{
Function (or character scalar for special case); the \code{imager} combine function to apply. This can be any of the ones in ?imager.combine. The default of NULL means it will use the \code{average} function. As a special case you can use \option{imager_func} = 'quantile' (with the speech marks) which uses \code{colQuantiles} in the \code{matrixStats} package.
}
  \item{na.rm}{
Logical; should NA be ignored in images when stacking (default TRUE)? Only applied to \option{imager_func} that has this as an option (otherwise ignored).
}
  \item{weights}{
Numeric vector; weights to use when stacking (must be the same length as the \option{image_list}). Only used when \option{imager_func} is \code{wsum}.
}
  \item{increasing}{
Logical; if TRUE, sort in increasing order (default TRUE). Only relevant for the \code{parsort}, \code{parorder}, \option{code} \option{imager_func} \code{imager} combining functions.
}
  \item{probs}{
Numeric vector; quantile probabilites between 0-1. This is only used in the specific case of \option{imager_func} = 'quantile' (with the speech marks). Behind the scenes it uses \code{colQuantiles} in the \code{matrixStats} package.
}
  \item{ondisk}{
Logical; should stacking be done using on-disk data? In this case all objects in \option{image_list} should be an Rfits_pointer. This can be done easily using \code{Rfits_make_list}.
}
  \item{cores}{
Integer scalar; number of cores used.
}
  \item{multitype}{
Character scalar; options are 'fork' (often best for Unix like and non-Gui sessions) and 'cluster' (Windows, and when using a Gui like R-Studio or threaded task). 'cluster' is also meant to be a bit more stable (if you have issues) but there is more memory copying overhead, so will generally be slower.
}
  \item{chunk}{
Integer scalar; chunk sizes to maximally break the image up into to compute the median stacks. The sweet spot seems to be around 1000 - 2000 depending on sytem RAM and number of \option{cores}.
}
}
\details{
The \code{propaneStackFlatInVar} stacking is actually done based on inverse variance weighting. In pseudo code:

stack = 0
stackRMS = 0
for(i in 1:length(image_list)){
  stack = stack + (image_list[[i]] - sky_list[[i]])/(skyRMS_list[[i]]^2)
  sky_stack = sky_stack + (image_list[[i]]^2)
}
stack = stack*sky_stack/(length(skyRMS_list)^2)

The output is explictly sky subtracted (so the sky is now 0 everywhere by definition as far as \code{profoundProFound} is concerned). The stacked sky is not returned. However, it can be computed by running \code{propaneStackFlatInVar} again, but passing the sky list originally passed to the \option{sky_list} argument to the \option{image_list} argument instead, and not providing any input to the \option{sky_list} argument (or setting this to 0).
}
\value{
For \code{propaneStackFlatInVar} a list containing:
  \item{image}{Numeric matrix; the variance-weighted sky-subtracted stacked image. Masked pixels are NA.}
  \item{skyRMS}{Numeric matrix/scalar; the sky RMS image/value of the final stacked image}
  \item{magzero}{The mag-zero point of the stacked image.}

For \code{propaneStackFlatFunc} a list containing:
  \item{image}{Numeric matrix; the combined image. Masked pixels are NA.}
  \item{weight}{Integer matrix; the pixel weights (number of non-masked pixels considered). Only relevant if \option{na.rm} is TRUE, and \code{imager} function being used has that as an option.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{propaneStackWarp}}
}
\examples{
library(ProFound)
library(Rfits)
library(magicaxis)

im1 = im2 = im3 = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))$imDat

stack = propaneStackFlatInVar(list(im1, im2, im3),
skyRMS_list = list(8,8,3))

#The new signal-to-noise weighted sky should equal sqrt(1/(1/8^2+1/8^2+1/3^2)) = 2.65

stack$skyRMS

# masking logic, here we have a wedding cake of masked regions:

im1[100:200,100:200]=NA; im2[120:180,120:180]=NA; im3[140:160,140:160]=NA

# masking='or' will conservatively mask any pixel that is masked in the stack:

magimage(propaneStackFlatInVar(list(im1, im2, im3), masking='or')$image)

# masking='and' will optimistically only mask pixels masked in all stacked images:

magimage(propaneStackFlatInVar(list(im1, im2, im3), masking='and')$image)

# this is the default manner for stacking in the more general propaneStackFlatFunc:

magimage(propaneStackFlatFunc(list(im1, im2, im3))$image)
magimage(propaneStackFlatFunc(list(im1, im2, im3))$weight, magmap=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ ~stack }% use one of  RShowDoc("KEYWORDS")
