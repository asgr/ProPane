\name{propaneInterp2D}
\alias{propaneInterp2D}
\alias{propaneBin2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rapidly Interpolate and Bin onto an Image
}
\description{
Rcpp function to rapidly interpolate or bin scatter data onto a target \option{image}. This is much faster than \code{akima} and similar functions (but probably less accurate). It at least guarantees flux is conserved (but maybe slightly in the wrong place).
}
\usage{
propaneInterp2D(x, y, z = 1, image, xlim = NULL, ylim = NULL, pixcen='FITS', type = 'add',
  zero = FALSE)
propaneBin2D(x, y, z = 1, image, xlim = NULL, ylim = NULL, pixcen='FITS', type = 'add',
  zero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector; x locations in FITS/R coordinates. See \option{pixcen}
}
  \item{y}{
Numeric vector; y locations in FITS/R coordinates. See \option{pixcen}
}
  \item{z}{
Numeric vector; z values to interpolate and add/subtract onto output \option{image} based on \option{type}. If length 1 then the same value is used for all \option{x} and \option{y} locations.
}
  \item{image}{
Numeric matrix; image to interpolate onto. Interpolated values will be added or subtracted from this based on \option{type}.
}
  \item{xlim}{
Numeric vector; x location of the first and last pixel in \option{image} in the same scheme as \option{pixcen}.
}
  \item{ylim}{
Numeric vector; y location of the first and last pixel in \option{image} in the same scheme as \option{pixcen}.
}
  \item{pixcen}{
Character scalar; either 'FITS' (pixel centres are integer, e.g. the centre of the bottom-left pixel is [1,1]) or 'R' (pixel centres are half integer, e.g. the centre of the bottom-left pixel is [0.5,0.5]). The default is 'FITS', to be consistent with other programs, but within R images are usually interpretted with a half pixel offset, so for converting object image [x,y] positions to true WCS positions option 'R' should be used.
}
  \item{type}{
Character scalar; if type = 'add' then interpolated values are added to the provided \option{image}, if type = 'sub' then interpolated values are subtracted from the provided \option{image}
}
  \item{zero}{
Logical; should the provided \option{image} be zeroed before adding/subtracting?
}
}
\details{
\code{propaneInterp2D} uses a simple linear weighting scheme of the fractional distance to the pixel boundary. Note this system assumes the data being intepolated inherently has the same orientation as the target pixel grid and pixels the same size.

\code{propaneBin2D} entirely places the data into the closest bin. This means it can operate as a rapid 2D histogram method.
}
\value{
Numeric matrix; the interpolated \option{image} data if not \option{xlim} and/or \option{ylim}. Otherwise a list with components: x (vector of x grid locations); y (vector of y grid locations); z (numeric matrix).
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{propaneWarp}}
}
\examples{
library(magicaxis)

magimage(propaneInterp2D(1:10, 1:10, 1:10, matrix(0,10,10)), qdiff=TRUE)
magimage(propaneInterp2D(1:10 - 0.25, 1:10 + 0.35, 1:10, matrix(0,10,10)), qdiff=TRUE)
magimage(propaneInterp2D(1:10 - 0.25, 1:10 + 0.95, -4:5, matrix(0,10,10)), qdiff=TRUE)
magimage(propaneInterp2D(runif(1e4,0,80), runif(1e4,0,80), 1, matrix(0,100,100)))

#We can also map the points onto an arbitrary grid by giving x and y limits:
magimage(propaneInterp2D(runif(1e3,-10,10), runif(1e3,-10,10), 1, image=matrix(0,100,100),
  xlim=c(-12.5,12.5), ylim=c(-12.5,12.5), pixcen='R'))

#Notice we conserve flux whether interpolating or binning
temp1 = propaneInterp2D(runif(1e3,-10,10), runif(1e3,-10,10), 1, image=matrix(0,100,100),
  xlim=c(-10,50), ylim=c(-10,50), pixcen='R')
sum(temp1$z)
magimage(temp1)

temp2 = propaneBin2D(runif(1e3,-10,10), runif(1e3,-10,10), 1, image=matrix(0,100,100),
  xlim=c(-10,50), ylim=c(-10,50), pixcen='R')
sum(temp2$z)
magimage(temp2)
}

\concept{ ~interpolation }

