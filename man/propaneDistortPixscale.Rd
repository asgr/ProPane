\name{propaneDistort}
\alias{propaneDistortPixscale}
\alias{propaneDistortPixarea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pixel Distortion Maps
}
\description{
Calculate a pixel distortion map for a full image, either in terms of the local pixel scale, or the local pixel area. This can be useful when consider flux and/or surface brightness in very large images where distortion can be non-negligible in the corners.
}
\usage{
propaneDistortPixscale(keyvalues, useraw = TRUE, unit = "asec", ...)
propaneDistortPixarea(keyvalues, useraw = TRUE, unit = "asec2", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keyvalues}{
List; standard \code{Rfits} list format key values header for the WCS of interest.
}
  \item{useraw}{
Logical; should the full raw header be used and passed to CFITSIO? This is often good (uses all the distortion terms etc), but there are failing cases with long messy headers so sometimes this cannot be used.
}
  \item{unit}{
Character scalar; for \code{propaneDistortPixscale} must be one of asec / amin / deg / rad, for \code{propaneDistortPixarea} must be one of asec2 / amin2 / deg2 / rad2 / str (where rad2 and str are both steradians).
}
  \item{\dots}{
Other arguments to pass to \code{Rwcs_p2s}.
}
}
\details{
Both functions return a full image, where pixel values are either the exact pixel scale (\code{propaneDistortPixscale}) or area (\code{propaneDistortPixarea}) of the relevant pixel. These full matrix outputs can be relatively expensive to compute for very big images, so consider using the \code{Rfits::pixscale} and \code{Rfits::pixarea} \code{Rfits} methods where you can inspect the extreme corners with the 'loc' option.
}
\value{
An object of class Rfits_image, where the imDat part is either the pixel scale or pixel area as requested. This will match 1-1 with the original image.
}
\author{
Aaron Robotham
}

\seealso{
\code{Rfits::pixscale}, \code{pixarea}
}
\examples{
\dontrun{
library(Rfits)
library(Rwcs)

VST_r = Rfits_read_image(system.file('extdata/VST_r.fits', package="magicaxis"))

VST_r_pixscale = propaneDistortPixscale(VST_r$keyvalues)

#check

VST_r_pixscale$imDat[1,1] == pixscale(VST_r, loc='bl')

extreme_pixscale = range(VST_r_pixscale$imDat)
extreme_pixscale[2] / extreme_pixscale[1] #This is the relative extreme variation
plot(VST_r_pixscale)
}
}
