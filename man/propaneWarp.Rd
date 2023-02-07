\name{propaneWarp}
\alias{propaneWarp}
\alias{propaneRebin}
\alias{propaneTran}
\alias{propaneWarpProPane}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remap and Rebin Image WCS via Warping
}
\description{
Remaps an input projection system to a different target WCS.

\code{propaneWarp} does precise remapping of one WCS to another, \code{propaneRebin} does quick coarse down-sampling (via re-binning), which can be useful for faster visualisations (since fewer pixels to worry about). Both conserve flux correctly. \code{propaneTran} just applies a translational shift. \code{propaneWarpProPane} warps all elements of a ProPane class object (as returned by e.g. \code{\link{propaneStackWarpInVar}}).
}
\usage{
propaneWarp(image_in, keyvalues_out = NULL, keyvalues_in = NULL, dim_out = NULL,
  pixscale_out = NULL, pixscale_in = NULL, direction = "auto", boundary = "dirichlet",
  interpolation = "cubic", doscale = TRUE, dofinenorm = TRUE, plot = FALSE,
  header_out = NULL, header_in = NULL, dotightcrop = TRUE, keepcrop = FALSE,
  WCSref_out = NULL, WCSref_in = NULL, magzero_out = NULL, magzero_in = NULL,
  blank = NA, warpfield = NULL, warpfield_return = FALSE, ...)

propaneWarpProPane(propane_in, keyvalues_out = NULL, dim_out = NULL,
  magzero_out = NULL, ...)

propaneRebin(image, scale = 1, interpolation = 6)

propaneTran(image, delta_x = 0, delta_y = 0, direction = 'backward', padNA = TRUE,
  shift_int = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_in}{
Numeric matrix; required, the image we want to warp. If \option{image} is a list as created by \code{Rfits_read_image} then the image part of the list is parsed to \option{image} and the correct header part is parsed to \option{keyvalues}.
}
  \item{propane_in}{
List; object of class 'ProPane', e.g. the output of \code{\link{propaneStackWarpInVar}}. All present list elements will be warp as appropriate (so conserving flux, inverse-Variance, or value as required): image (flux), weight (value), inVar (rebinned inVar), exp (value), cold (flux), hot (flux), clip (value).
}
  \item{image}{
Numeric matrix; required, the image we want to rebin. If \option{image} is a list as created by \code{Rfits_read_image} then the image part of the list is parsed to \option{image} and the correct header part is parsed to \option{keyvalues}.
}
  \item{keyvalues_out}{
List; output header values to be used for the WCS. This is the target WCS projection that \option{image_in} will be mapped onto.
}
  \item{keyvalues_in}{
List; input header values to be used for the WCS. This should be the header WCS that matches \option{image_in}.
}
  \item{dim_out}{
Integer vector; this defines the desired dimensions of the output image. If this is not provided then the output image is made to be the same size as the NAXIS1 and NAXIS2 arguments taken from \option{header_out} (which is usually what you will want TBH).
}
  \item{pixscale_out}{
Numeric scalar; output pixel scale in asec/pixel. If not provided it will be computed from the \option{keyvalues_out} using \code{\link{Rwcs_pixscale}}, which is only an issue if the CD matrix is defined a long way from the actual projected image which can lead to distorted pixel scales.
}
  \item{pixscale_in}{
Numeric scalar; input pixel scale in asec/pixel. If not provided it will be computed from the \option{keyvalues_in} using \code{\link{Rwcs_pixscale}}, which is only an issue if the CD matrix is defined a long way from the actual projected image which can lead to distorted pixel scales.
}
  \item{direction}{
Character scalar; "auto" (default), "forward" or "backward", see \code{imwarp}. Since it is usally better to go from the higher resolution image and map this onto the lower resolution grid, "auto" selects the better direction given the pixel scales recovered from the header information.
}
  \item{boundary}{
boundary conditions: "dirichlet", "neumann", "periodic" (default "dirichlet"), see \code{imwarp}
}
  \item{interpolation}{
Character/integer scalar; For \code{propaneWarp}: "nearest", "linear", "cubic" (default "linear"), see \code{imwarp}. For \code{imresize}: Linear is 3, Cubic is 5, Lanczos (higher quality), default) is 6.
}
  \item{doscale}{
Logical; if TRUE (default) then the image is scaled by the change in pixel scale, i.e. it properly conserves angular flux density. This is usually what you minimally want when warping a flux image between different WCS. If FALSE then nothing is done to rescale the image. This is useful when the input \option{image_in} is something like a segmentation or exposure map, where the actual values matter.
}
  \item{dofinenorm}{
Logical; if TRUE (default) then the image is scaled by the relative change the change in flux per pixel found in a reference normalisation image. This is usually what you want when converting images between different WCS (since we really want to make sure we conserve overall flux). It mostly adjusts fluxes around the image edge (where interpolation artefacts may cause flux loss) and changes things often by on milli-mag. If that level of accuracy is not required then turning this to FALSE will speed things up considerably (by almost a factor of 2 for big images).
}
  \item{plot}{
Logical; should a \code{\link{Rwcs_image}} plot of the output be generated?
}
  \item{header_out}{
Optional output header. To aid backwards compatibility, you can also pass in a 1D character vector of header values, e.g. the \option{hdr} component of \option{readFITS} from the \code{FITSio} package. \option{keyvalues_out} must be left as NULL in this case. If you want to use all distortion terms you can also pass in the raw header output of \code{Rfits_read_header_raw}.
}
  \item{header_in}{
Optional input header. To aid backwards compatibility, you can also pass in a 1D character vector of header values, e.g. the \option{hdr} component of \option{readFITS} from the \code{FITSio} package. \option{keyvalues_in} must be left as NULL in this case. If you want to use all distortion terms you can also pass in the raw header output of \code{Rfits_read_header_raw}.
}
  \item{dotightcrop}{
Logical; should an internal tight crop be made. If the images are approximately the same angular size this should be left as FALSE (it will run faster), but when there is a big difference between the input and output WCS (the latter being much larger) then setting this to TRUE might speed the warping up a lot.
}
  \item{keepcrop}{
Logical; if \option{dotightcrop} is TRUE, then this flag specifies if just the minimally cut down projection should be returned (TRUE; using less memory, but perhaps harder to work with) or the output should be fully cast onto the target WCS (FALSE).
}
  \item{WCSref_out}{
Scalar; optional output \option{WCSref}. Either NULL for base WCS (default), or a letter [a-z] for the relevant projection desired. Only used if using a raw \option{header}, and obviously this has to exist within the provided raw \option{header}.
}
  \item{WCSref_in}{
Scalar; optional input \option{WCSref}. Either NULL for base WCS (default), or a letter [a-z] for the relevant projection desired. Only used if using a raw \option{header}, and obviously this has to exist within the provided raw \option{header}.
}
  \item{magzero_out}{
Numeric scalar; the output mag-zero point desired. Note \option{magzero_in} must also be provided for any scaling to be applied. A nice \option{magzero_out} is perhaps 23.9, because this means all measured fluxed are in micro-Jansky (rather than arbitrary ADUs). \option{magzero_out} = 8.9 gives Jansky; 16.4 gives milli-Jansky and 31.4 gives nano-Jansky. For \code{propaneWarpProPane}, leaving as NULL means it is inherited from the MAGZERO keyword attached to the \option{promo_in} header.
}
  \item{magzero_in}{
Numeric scalar; the input mag-zero point. Note \option{magzero_out} must also be provided for any scaling to be applied.
}
  \item{scale}{
Numeric scalar; resolution increasing factor, i.e. value larger than 1 will increase the number of output pixels (decrease the asec/pix pixel scale) and vica-versa.
}
  \item{delta_x}{
Numeric scalar; the shift in apply in the x-axis of the input \option{image}. Positive means a visual shift to the right (with padding applied to the left).
}
  \item{delta_y}{
Numeric scalar; the shift in apply in the y-axis of the input \option{image}. Positive means a visual shift up (with padding applied to the bottom).
}
  \item{padNA}{
Logical; should vacant pixels be padded with NA (TRUE) or left as 0 (FALSE).
}
  \item{shift_int}{
Logical; should the shift applied be purely integer (TRUE) or fractional pixel (FALSE). Integer shifts are much faster to compute, and are better at conserving flux and not introducing pixel covariance, so think carefully before setting to FALSE.
}
  \item{blank}{
Scalar; the value to initially fill the output blank matrix with. In effect, this means that regions that receive to warped data (empty/blank regions) will all inherit this value. Sensibly you might want this to be NA (the default, implying missing), or for stacking purposes 0.
}
  \item{warpfield}{
Cimg imager_array; a warp field to be used. You almost very certainly never want to actually make this yourself, but rather pass in a previously saved \option{warpfield} output for speeding up future warping operations.
}
  \item{warpfield_return}{
Logical; should a \option{warpfield} be appended to the output list? This can be used later with the \option{warpfield} argument to speed up identical WCS warping operations.
}
  \item{\dots}{
For \code{propaneWarp} dots are passed to \code{\link{Rwcs_image}} (only relevant if \option{plot}=TRUE). For \code{propaneWarpProPane} dots are passed to \code{propaneWarp}.
}
}
\details{
The \code{propaneWarp} function allows for arbitrary WCS remapping from any \code{wcslib} supported projection (see \code{\link{Rwcs}}). The process internally does the following:

\itemize{
  \item Rwcs_p2s; maps the input image to RA and Dec per pixel using the input header
  \item Rwcs_s2p; maps the pixel RA and Dec coordinates onto the desired output pixel grid using the output header
}

There are a few different ways to consider doing pixel remapping (or warping). The main question is usually whether to operate on a forward or backwards manner (see \code{imwarp}). Backwards mapping finds the best interpolation of every output pixel in the input image, whereas forward mapping finds the best interpolation of every input pixel on the output image. Backwards mapping usually provides the fewest artefacts, but in general if the input image is higher resolution then you might prefer forward mapping and vica-versa.

The actual warping is done using the \code{imwarp} function in the \code{imager} package, so this will need to be installed prior to use (it is available on CRAN).

The \code{propaneRebin} allows simple rebinning on the same projection geometry (so no rotation or distortion term differences). For this reason it is simpler to use, and also work on matrix inputs (really with PSF re-scaling in mind).
}
\value{
For both \code{propaneWarp} and \code{propaneRebin} (if input is class 'Rfits'), returns a list (class 'Rfits_image') containing:

  \item{imDat}{Numeric matrix; the remapped image using the target WCS.}
  \item{keyvalues}{Usual keyvalues of target \option{keyvalues_out}.}
  \item{hdr}{Usual hdr of target \option{header_out}.}
  \item{header}{Usual header of target \option{header_out}.}
  \item{header}{Usual raw of target \option{header_out}.}
  \item{keynames}{Usual keynames of target \option{header_out}.}
  \item{keycomments}{Usual keycomments of target \option{header_out}.}

For \code{propaneRebin}, if the input is a matrix then the output is also a matrix. The user will need to pay attention to the change in the pixel scale etc because this is therefore not stored directly in the output since there is no attached WCS.

For \code{propaneWarp} if \option{dotightcrop} and \option{keepcrop} are both TRUE then the output will be tightly cropped and the limits of the image subset this represents in the target output image exist in an additional \option{crop} list output (a vector of [xlo,xhi,ylo,yhi]). This can be used to insert the projected image in the correct place in a larger target image (potentially on disk). Regardless, the WCS relative to the cropped pixels will still be correct.

\code{propaneTran} simply returns the shifted numeric matrix.

\code{propaneWarpProPane}

If \option{warpfield_return} is TRUE then \code{propaneWarp} will also have a \option{warpfield} list output. The idea is this can be passed back into the \option{warpfield} argument of a future \code{propaneWarp} operation to speed up computation.

The output is an object of class 'ProPane' (warped as appropriate) containing:

\item{image}{The stacked image, of class Rfits_image.}
\item{weight}{The number of available un-masked pixels in every output pixel of \option{image}. Class Rfits_image.}
\item{inVar}{If \option{inVar_list} is not NULL, the stacked inverse variance of the output \option{image}. Class Rfits_image when present.}
\item{exp}{If \option{exp_list} is not NULL, the total exposure time of the output \option{image}. Class Rfits_image when present. Note this will be unmasked, so values of 0 means no photon ever hit the output \option{image}, i.e. masked pixels can have postive exposures if they sit within the detector domain.}
\item{cold}{If \option{keep_extreme_pix} = TRUE, the extreme hot pixel value for the input \option{image_list} that formed the final stacked \option{image}. If \option{keep_extreme_pix} = FALSE this will be NULL. Class Rfits_image when present.}
\item{hot}{If \option{keep_extreme_pix} = TRUE, the extreme hot pixel value for the input \option{image_list} that formed the final stacked \option{image}. If \option{keep_extreme_pix} = FALSE this will be NULL. Class Rfits_image when present.}
\item{clip}{The number of frames clipped from the stacked \option{image}. This is only relevant when \option{doclip} is TRUE. Values larger than 1 can indicate that both cold and hot pixels have been clipped in this location, and also that the dilation process has masked a pixel multiple times.}

}
\author{
Aaron Robotham
}
\note{
This function uses a bi-cubic interpolation scheme by default. It should approximately conserve the flux in \option{image_in}, but this is not guaranteed. The 'linear' interpolation scheme may be closer to conserving flux in images with sharp features. The conservation is usually good to about 0.1\% (i.e. 0.01 mag). If you require better conservation than this then alternative schemes (e.g. SWarp) will be required.

For working with segmentation maps, the following options should usually be used: \option{doscale} = FALSE, \option{interpolation} = 'nearest', \option{direction} = 'backward'.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Rwcs}}, \code{\link{propaneStackWarp}}
}
\examples{
\dontrun{
library(Rfits)
GALEX_NUV=Rfits_read_image(system.file("extdata", 'GALEX_NUV.fits', package="magicaxis"))
VST_r=Rfits_read_image(system.file("extdata", 'VST_r.fits', package="magicaxis"))
VISTA_K=Rfits_read_image(system.file("extdata", 'VISTA_K.fits', package="magicaxis"))

propaneWarp(VST_r, GALEX_NUV$keyvalues, plot=TRUE)
propaneWarp(VISTA_K, GALEX_NUV$keyvalues, plot=TRUE)

propaneWarp(GALEX_NUV, VST_r$keyvalues, plot=TRUE)
propaneWarp(VISTA_K, VST_r$keyvalues, plot=TRUE)

propaneWarp(GALEX_NUV, VISTA_K$keyvalues, plot=TRUE)
propaneWarp(VST_r, VISTA_K$keyvalues, plot=TRUE)

#Check we can warp forwards and backwards correctly:

plot(GALEX_NUV)
propaneWarp(GALEX_NUV, GALEX_NUV$keyvalues, plot=TRUE)
propaneWarp(propaneWarp(GALEX_NUV, VST_r$keyvalues, interpolation='nearest'),
        GALEX_NUV$keyvalues, interpolation='nearest', plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ ~warp }% use one of  RShowDoc("KEYWORDS")