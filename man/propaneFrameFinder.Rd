\name{propaneFrameFinder}
\alias{propaneFrameFinder}
\alias{propaneGenWCS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finder Frame Coverage
}
\description{
For a given region of sky and a number of on-disk FITS files, the \code{propaneFrameFinder} function will find all nearby frames and optional plot their coverage. The \code{propaneGenWCS} can compute the most efficient WCS that covers the target frames. This is useful when trying to figure out an optimal WCS for stacking purposes.
}
\usage{
propaneFrameFinder(filelist = NULL, dirlist = NULL, extlist = 1, RAcen = 180,
  Deccen = 0, rad = 1, cores = 1, plot = TRUE, proj = TRUE, border = "red",
  col = hsv(alpha = 0.2), bg = par("bg"), ...)

propaneGenWCS(filelist = NULL, dirlist = NULL, image_list = NULL, rotation = 'North',
  CRVAL1 = NULL, CRVAL2 = NULL, pixscale = NULL, NAXIS1 = NULL, NAXIS2 = NULL,
  CRPIX1 = NULL, CRPIX2 = NULL, CTYPE1 = "RA---TAN", CTYPE2 = "DEC--TAN",
  CD1_1 = NULL, CD1_2 = NULL, CD2_1 = NULL, CD2_2 = NULL, CUNIT1 = "deg",
  CUNIT2 = "deg", ...)
}
\arguments{
  \item{filelist}{
Character vector; vector of full paths of FITS files to analyse. Both \option{filelist} and \option{dirlist} can be provided, and the unique superset of both is scanned.
}
  \item{dirlist}{
Character vector; vector of full paths of directories to analyse. Both \option{filelist} and \option{dirlist} can be provided, and the unique superset of both is scanned.
}
  \item{extlist}{
Integer vector; the extensions to use. If length 1 then it will be used for all \option{filelist}, but otherwise it should be the same length as the final \option{filelist}
}
  \item{RAcen}{
Numeric scalar; the Right Ascension of the centre of the field in degrees.
}
  \item{Deccen}{
Numeric scalar; the Declination of the centre of the field in degrees.
}
  \item{rad}{
Numeric scalar; the radius from \option{RAcen}, \option{Deccen} to search in units of degrees.
}
  \item{cores}{
Integer scalar; the number of cores to run on.
}
  \item{plot}{
Logical; should a coverage plot be generated?
}
  \item{proj}{
If \option{plot} = TRUE then should \code{magproj} be used for more slower but more accurate projection?
}
  \item{border}{
Scalar; colour of the border of frames when \option{plot} = TRUE.
}
  \item{col}{
Scalar; colour inside of frames when \option{plot} = TRUE.
}
  \item{bg}{
Scalar; the background color for the legend box (note on some devices [e.g. Quartz] the default it solid white, and others [e.g. pdf] transparent).
}
  \item{image_list}{
List; required. List of the input \code{Rfits_image} or \code{Rfits_pointer} objects to be warped and stacked. If this is not present then one or both of \option{filelist} and/or \option{dirlist} must be provided.
}
  \item{rotation}{
Character/numeric scalar; the alignment of the image. If character then can be the direction which is up, i.e. 'North', 'East', 'South', 'West', or 'get'. If 'get' then it will determine the average rotation of all the target images. If numeric then it specifies the manual rotation to apply. This is the vector of North with respect to y-axis up clockwise (so -> is 90, i.e. East is up).
}
  \item{CRVAL1}{
Numeric scalar; FITS header CRVAL1 for the CTYPE1 projection system. This is the RA in degrees at the location of \option{NAXIS1}/2. This should be chosen to be near the centre of the target projection (i.e. the average RA of the frames to be overlaid). If not provided this will be estimated from the provided data.
}
  \item{CRVAL2}{
Numeric scalar; FITS header CRVAL2 for the CTYPE2 projection system. This is the Dec in degrees at the location of \option{NAXIS2}/2. This should be chosen to be near the centre of the target projection (i.e. the average Dec of the frames to be overlaid). If not provided this will be estimated from the provided data.
}
  \item{pixscale}{
Numeric scalar; the target pixel scale in arc seconds. By default \option{CD1_1} and \option{CD2_2} inherit this in degress (i.e. divided by 3600). Note the standard convention is to have RA increasing as you move left on an image, therefore \option{CD1_1} = -\option{pixscale}/3600 by default (as seen in the arguments). If not provided this will be estimated from the provided data.
}
  \item{NAXIS1}{
Numeric scalar; the number of pixels to display on the x-axis. If not provided this will be estimated from the provided data.
}
  \item{NAXIS2}{
Numeric scalar; the number of pixels to display on the y-axis. If not provided this will be estimated from the provided data.
}
  \item{CRPIX1}{
Numeric scalar; FITS header CRPIX1 for the CTYPE1 projection system. This is the x pixel value at the location of CRVAL1. The default is NAXIS1/2 + 0.5, i.e. the centre of the image (with the FITS definition of centres of pixels being integer values). If not provided this will be estimated from the provided data.
}
  \item{CRPIX2}{
Numeric scalar; FITS header CRPIX2 for the CTYPE2 projection system. This is the y pixel value at the location of CRVAL2. The default is NAXIS2/2 + 0.5, i.e. the centre of the image (with the FITS definition of centres of pixels being integer values). If not provided this will be estimated from the provided data.
}
  \item{CTYPE1}{
Character scalar; the RA projection system type. Must be a WCS legal type (up to 4 letter axis followed by 3 letter projection padded as required with "-", e.g.: "RA---TAN") and 8 characters long. If not provided this will be estimated from the provided data.
}
  \item{CTYPE2}{
Character scalar; the Dec projection system type. Must be a WCS legal type (up to 4 letter axis followed by 3 letter projection padded as required with "-", e.g.: "DEC--TAN") and 8 characters long. If not provided this will be estimated from the provided data.
}
  \item{CD1_1}{
Numeric scalar; FITS header CD1_1 for the CTYPE1 projection system. Change in CTYPE1 in degrees along x-Axis. If not provided this will be estimated from the provided data.
}
  \item{CD1_2}{
Numeric scalar; FITS header CD1_2 for the CTYPE1 projection system. Change in CTYPE1 in degrees along y-Axis. If not provided this will be estimated from the provided data.
}
  \item{CD2_1}{
Numeric scalar; FITS header CD2_1 for the CTYPE2 projection system. Change in CTYPE2 in degrees along x-Axis. If not provided this will be estimated from the provided data.
}
  \item{CD2_2}{
Numeric scalar; FITS header CD2_2 for the CTYPE2 projection system. Change in CTYPE2 in degrees along y-Axis. If not provided this will be estimated from the provided data.
}
  \item{CUNIT1}{
Character scalar; the units of the WCS for the RA axis. Is de facto 'deg' almost always these days, to the degree that many headers do not even list it (or list it and leave it blank) assuming the user will assume it is in degrees. If found and set to a blank string it will be changed to 'deg'. If not provided this will be estimated from the provided data.
}
  \item{CUNIT2}{
Character scalar; the units of the WCS for the Dec axis. Is de facto 'deg' almost always these days, to the degree that many headers do not even list it (or list it and leave it blank) assuming the user will assume it is in degrees. If found and set to a blank string it will be changed to 'deg'. If not provided this will be estimated from the provided data.
}
  \item{\dots}{
Other arguments to pass to \code{\link{Rfits_key_scan}}.
}
}
\details{
At the simplest level this function uses \code{\link{Rfits_key_scan}} to get useful information about target FITS files. It then uses \code{\link{coordmatchsing}} to find frame centres that are within the circle specified by \option{RAcen}, \option{Deccen} and \option{rad} and trims the output of \code{\link{Rfits_key_scan}} to just these frames. Optionally it can then plot the coverage.
}
\value{
For \code{propaneFrameFinder} the output is the trimmed output of \code{\link{Rfits_key_scan}} plus a separation column (sep) in units of \option{radunit}.

For \code{propaneGenWCS} the output is the automatic keyvalues list that best covers the provided data. Note this might slightly over-extend the minimal pixel limits needed (NAXIS1 and NAXIS2) by 1 or 2 pixels because the RA and Dec limits are estimated for a worst-case scenario projection.
}
\author{
Aaron Robotham
}
\note{
It is not necessary to trim down inputs before running \code{\link{propaneStackWarpInVar}}, but this function makes it easy to find an optimal WCS for the stacked frames.
}
\seealso{
\code{\link{Rfits_key_scan}}
}
\examples{
dir_stack = system.file('extdata/stack/', package="ProPane")

ex_frames = propaneFrameFinder(dirlist=dir_stack, extlist = 2, RAcen=36.9, Deccen=-5.2,
  rad=0.1, plot=TRUE)

#Note we will also find frames if our search region sits entirely inside them:

ex_frames = propaneFrameFinder(dirlist = dir_stack, extlist = 2, RAcen=36.85, Deccen=-5.2,
  rad=0.01, plot=TRUE)

print(ex_frames)

\dontrun{
image_list = Rfits_make_list(dirlist = dir_stack, extlist = 2, pointer = FALSE)

keyvalues_out = propaneGenWCS(image_list = image_list)

stack = propaneStackWarpInVar(image_list = image_list,
                   inVar_list = 0.000468, #assume all have the same inVar
                   exp_list = 10, #all have the same exposure time of 10s
                   cores = 8,
                   keyvalues_out = keyvalues_out,
                   magzero_in = 30,
                   magzero_out = 23.9 #micro-jansky output
                   )

plot(stack$image)
}
}

