\name{propaneCovTest}
\alias{propaneCovTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check Image Covariance
}
\description{
Allows users to experiment with the impact of different warping options by creating a sparse checkboard of pixels and computing how this discrete flux is re-redistributed.
}
\usage{
propaneCovTest(image_in, keyvalues_out = NULL, sparse = 5L, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_in}{
List; required, the image we want to warp. \option{image_in} must be either object class Rfits_image.
}
  \item{keyvalues_out}{
List; output header values to be used for the WCS. This is the target WCS projection that \option{image_in} will be mapped onto.
}
  \item{sparse}{
Integer scalar; how sparse we want the checkerboard pattern to be. Generally sparser is better so you can be certain sparse pixels do not interfer with each other, but you still want it to be dense enough to capture small scale features in the warping.
}
  \item{\dots}{
Other arguments passed on to \code{\link{propaneWarp}}.
}
}
\details{
The checkerboard grid is project with the arguments as specified, and then the output image is assessed for how cleanly projected the grid is. The brighter (1 - sparse^2) pixels are assumed to be those that correspond most directly to the input coarse pixels, and these should ideally average 1 (\option{hi_mean}). The surrounding 8 projected pixels are also assessed for local covariance (at the single pixel leve), and these should ideally average 0 (\option{lo_mean}).
}
\value{
\item{cov_im }{Rfits_image; the warped coarse pixel image.}
\item{lo_mean }{The mean value of the 8 projected pixels around each input coarse pixel (nearer to 0 is better).}
\item{hi_mean }{The mean value of the projected coarse pixel (nearer to 1 is better).}
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{propaneWarp}}
}

\examples{
\dontrun{
image_list = Rfits_make_list(dirlist = system.file('extdata/stack/', package="ProPane"),
                              extlist = 2, pointer = FALSE) #load full images

image_1 = image_list[[1]][,]
image_2 = image_list[[2]][,]

cov_check = propaneCovTest(image_2, image_1$keyvalues, direction='backward')

plot(cov_check$cov_im[1:100,1:100], mgamp=FALSE)

cov_check$lo_mean #This is 0.009, so pretty low (good)
cov_check$hi_mean #This is 0.936, so pretty high (good)
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
