\name{propaneLocalFunc}
\alias{propaneLocalFunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computer Statistics on Local Pixels
}
\description{
Computes local quantites such as medians but on a very small dithered regions around each pixel. The default computes statistics on a 3x3 local grid.
}
\usage{
propaneLocalFunc(image, imager_func = NULL, dither = 1, offset = 1, iter = 1,
  kern = 'square', cores = 1, multitype = "fork", verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; an image to process.
}
  \item{imager_func}{
Function (or character scalar for special case); the \code{imager} combine function to apply. This can be any of the ones in ?imager.combine. The default of NULL means it will use the \code{average} function. As a special case you can use \option{imager_func} = 'quantile' (with the speech marks) which uses \code{colQuantiles} in the \code{matrixStats} package. \code{parmed} will use the parallel median combine (which is one of the more useful).
}
  \item{dither}{
Integer scalar/vector; the number of dithers to execute, each one will be of size \option{sparse}. If vector must be length \option{iter}, in which case the i^th iteration will use \option{dither}[i] as the dither value.
}
  \item{offset}{
Integer scalar/vector; the number of pixels used for each dither. If vector must be length \option{iter}, in which case the i^th iteration will use \option{offset}[i] as the offset value.
}
  \item{iter}{
Integer scalar; the number of iterations to apply. Each new iteration will run \code{propaneLocalFunc} on the output of the previous iteration. If \option{dither} / \option{offset} is a vector The i^th iteration will use \option{dither}[i] / \option{offset}[i] as the dither / offset value.
}
  \item{kern}{
Character scalar; the type of local kernel to use. Options are 'square' which comptues a full grid of local dithers (so a square pattern) or 'circle' (which only users dithers that lie within the radius of \option{dither}). \option{dither} = 'circle' uses fewer dithers at larger radii for this reason, and will be a bit faster to run (about 30\%).
}
  \item{cores}{
Integer scalar; number of cores used to create the dithered offsets.
}
  \item{multitype}{
Character scalar; options are 'fork' (often best for Unix like and non-Gui sessions) and 'cluster' (Windows, and when using a Gui like R-Studio or threaded task). 'cluster' is also meant to be a bit more stable (if you have issues) but there is more memory copying overhead, so will generally be slower.
}
  \item{verbose}{
Logical; should a verbose output be written to screen (showing the iteration number)?
}
  \item{\dots}{
Other arguments passed to \code{\link{propaneStackFlatFunc}}.
}
}
\details{
The core routine will compute medians for local subsets of every pixel. This will produce an image with the same dimensions as the original but now processed using the specified \option{imager_func}. If iterations are being used (via the \option{iter} argument) then the function is re-run on each subsequent output. This allows for quite complex cleaning and pixel averaging, allowing different type of control over the structure to remove (or not).}
\value{
Numeric matrix with dimensions of the input \option{image}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{propaneStackFlat}}
}
\examples{
\dontrun{
library(Rfits)
library(magicaxis)
library(imager)

image = Rfits_read_image(system.file('extdata/VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))$imDat
magimage(image, qdiff=TRUE)

magimage(propaneLocalFunc(image), qdiff=TRUE)

#Using local median:
temp_local = propaneLocalFunc(image, imager_func=parmed, dither=2L, iter=5L)
magimage(temp_local, qdiff=TRUE)
magimage(image - temp_local, qdiff=TRUE)

#Using local median and increasing offset (for each iteration gets larger):
temp_local = propaneLocalFunc(image, imager_func=parmed, dither=2L, offset=1:5, iter=5L)
magimage(temp_local, qdiff=TRUE)
magimage(image - temp_local, qdiff=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
