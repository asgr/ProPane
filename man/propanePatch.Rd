\name{propanePatch}
\alias{propanePatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Path Inverse Variance Weighted Images with Median Combined Images
}
\description{
This function offers a few routes to optimally patch an inverse variance weighted image with a median combined one of the same region. The reason you might do this is inverse variance weighted stacking will have a minimum factor 1.253 less RMS than the same median combined data (and that's in the regime where each image is equal depth, if it varies then median combining gets relatively noisier).

However, it can be hard to fully clip pixels when stacking many frames because the chance of two frames producing hot/cold pixels in the same position increases. Using this function we aim to get the best of both worlds- the maximal depth of inverse variance weighted stacking and the better pixel clipping properties of a median combined image, i.e. we aim to patch the bad hot/cold inverse variance pixels with the preferred median combined pixels in the same position.
}
\usage{
propanePatch(image_inVar, image_med, diff_type = "scale", threshold = 5,
  scale_type = "quan", scale_val = 0.3, hot = TRUE, cold = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_inVar}{
Rfits_image or numeric matrix; required, the inverse variance weighted image data.
}
  \item{image_med}{
Rfits_image or numeric matrix; required, the median combined image data.
}
  \item{diff_type}{
Character scalar; the type of difference to consider when patching pixels. Options are 'scale' (the default) where we patch pixels according to diff = (\option{image_inVar} - \option{image_med})/val (val depends on \option{scale_type} and \option{scale_val}); 'rel' (default) where we patch pixels according to diff = (\option{image_inVar} - \option{image_med})/max(\option{image_med}, val) (val depends on \option{scale_type} and \option{scale_val}), i.e. we consider the relative difference in the regime where image flux dominates (usually inside of objects); 'diff' where we patch pixels according to diff = \option{image_inVar} - \option{image_med} (so just the raw difference).
}
  \item{threshold}{
Numeric scalar; given the \option{diff_type} we then threshold pixels using diff > \option{threshold}.
}
  \item{scale_type}{
Character scalar; the method of rescaling the image differences (the denominator 'val' discussed in \option{diff_type}). Options are 'quan' (the default) where val = quantile(abs(\option{image_med}), probs=\option{scale_val}, na.rm=TRUE); 'median' where val = median(abs(\option{image_med}), na.rm=TRUE); 'val' or 'num' where val = \option{scale_val}.
}
  \item{scale_val}{
Numeric scalar; argument to use in combination with \option{scale_type}.
}
  \item{hot}{
Should hot pixels be clipped?
}
  \item{cold}{
Should cold pixels be clipped?
}
}
\details{
When using this function it is often worth experimenting away from default settings a bit. You usually would not want to be patching in more than a couple percent of median combined pixels, and if the printed message suggests you are patching much more than this then certainly check the \option{patch} output to confirm you are not being too aggressive.

In general \option{diff_type} = 'scale' seems to work best on a range of images tested, but detectors with different characteristics might prefer 'rel' or even 'diff' options. Again, it is worth experimenting for different data types.

Note NA or NaN values in \option{image_inVar} are always patched with values from \option{image_med}.
}
\value{
List containing:
\item{image }{Rfits_image or numeric matrix; patched output image.}
\item{patch }{Rfits_image or numeric matrix; integer image showing the patched pixels (1L means patched due to NA, 2L means patched due to \option{threshold}).}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{propaneStackWarp}}
}
\examples{
# Nothing yet!
}

