\name{propaneStackWarp}
\alias{propaneStackWarp}
\alias{propaneStackWarpInVar}
\alias{propaneWarpDump}
\alias{propaneStackWarpMed}
\alias{propaneStackWarpFunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Warping and Stacking
}
\description{
Stacks multiple un-aligned images. This required warping to a common WCS before doing a traditional flat stack.
}
\usage{
propaneStackWarpInVar(image_list = NULL, inVar_list = NULL, exp_list = NULL,
  weight_list = NULL, mask_list = NULL, magzero_in = 0, magzero_out = 23.9,
  keyvalues_out = NULL, dim_out = NULL, cores = floor(detectCores()/2), cores_warp = 1,
  multitype = 'fork', Nbatch = cores, keepcrop = TRUE, keep_extreme_pix = FALSE,
  doclip = FALSE, clip_tol = 100, clip_dilate = 0, clip_sigma = 5, return_all = FALSE,
  dump_frames = FALSE, dump_dir = tempdir(), ...)

propaneWarpDump(image_list = NULL, magzero_in = 0, magzero_out = 23.9,
  keyvalues_out = NULL, dim_out = NULL, cores = floor(detectCores()/2), cores_warp = 1,
  multitype = 'fork', keepcrop = TRUE, dump_dir = tempdir(), ...)

propaneStackWarpMed(filelist = NULL, dirlist = NULL, extlist = 1,
  pattern = NULL, recursive = TRUE, zap = NULL, keyvalues_out = NULL,
  cores = floor(detectCores()/2), multitype = 'fork', chunk = 1e3, doweight = TRUE,
  useCUTLO = TRUE)

propaneStackWarpFunc(filelist = NULL, dirlist = NULL, extlist = 1,
  pattern = NULL, recursive = TRUE, zap = NULL, keyvalues_out = NULL, imager_func = NULL,
  weights = NULL, prob = 0.5, cores = floor(detectCores()/2),
  multitype = 'fork', chunk = 1e3, doweight = TRUE, useCUTLO = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_list}{
List; required. List of the input \code{Rfits_image} or \code{Rfits_pointer} objects to be warped and stacked. Masked pixels can be flagged with NA in the imDat component. The data must already be background subtracted to the level desired by the user (this is not done internally). Note that the header information attached to the \option{image_list} is used for \option{inVar_list}, \option{exp_list}, \option{weight_list} and \option{mask_list} as appropriate. See Details below.
}
  \item{inVar_list}{
Numeric vector or list of numeric matrices / \code{Rfits_image} / \code{Rfits_pointer}; optional. The inverse-variances. Either one value per \option{image_list} entry, or one pixel matched matrix per \option{image_list} entry. This is used to achieve S/N weighted stacking. If not provided then images are stacked with equal weight. See Details below.
}
  \item{exp_list}{
Numeric vector or list of numeric matrices / \code{Rfits_image} / \code{Rfits_pointer}; optional. The exposure time (in whatever units you like, but they should be consistent between input images). Either one value per \option{image_list} entry, or one pixel matched matrix per \option{image_list} entry. This is used to track maximal pixel exposure (so nothing is masked). This is mostly provided to persist knowledge from a previous stack that you are now adding to with more data. See Details below.
}
  \item{weight_list}{
Integer vector or list of integer matrices / \code{Rfits_image} / \code{Rfits_pointer}; optional. The weight map showing how many images have contributed to a given pixel. Either one value per \option{image_list} entry, or one pixel matched matrix per \option{image_list} entry. This is mostly provided to persist knowledge from a previous stack that you are now adding to with more data (else the previous stack would be implictly given the weight of 1 everywhere). See Details below.
}
  \item{mask_list}{
List of integer matrices / \code{Rfits_image} / \code{Rfits_pointer}; optional. Pixel matched mask matrix per \option{image_list} entry. Anything non-zero is considered masked. This can be used in conjuction with setting pixels in \option{image_list} to NA. See Details below.
}
  \item{magzero_in}{
Numeric vector; the input mag-zero points. If length 1 then it is assumed all input frames have the same mag-zero point.
}
  \item{magzero_out}{
Numeric scalar; the output mag-zero point desired. A nice \option{magzero_out} is perhaps 23.9 (the default), because this means all measured fluxed are in micro-Jansky (rather than arbitrary ADUs). \option{magzero_out} = 8.9 gives Jansky; 16.4 gives milli-Jansky and 31.4 gives nano-Jansky.
}
  \item{keyvalues_out}{
List; output header values to be used for the WCS. This is the target WCS projection that each image in \option{image_list} will be mapped onto. For help creating this have a look at \code{\link{propaneGenWCS}}.
}
  \item{imager_func}{
Function (or character scalar for special case); the \code{imager} combine function to apply. This can be any of the ones in ?imager.combine. The default of NULL means it will use the \code{average} function. As a special case you can use \option{imager_func} = 'quantile' / 'quan' (with the speech marks) for quantiles (defined by \option{probs}). \option{imager_func} = 'InVar' (with speech marks) will return the inverse-variance scaled correctly by exposures (which is what we use for general stacking). \option{imager_func} = 'waverage' will return a weighted average using the \option{weights} values for the per frame weights, but note this only works in regions where all frames contribute. In most practical applications you should use \code{propaneStackWarpInVar} for proper inverse-variance weighted stacking, but this method might be useful as a quick check.
}
  \item{weights}{
Numeric vector; weights to use for image combining. Only relevant for \option{imager_func} = imager::wsum or \option{imager_func} = 'waverage'. The latter becomes a weighted average, where the input weights should in effect be the per image inverse-variance.
}
  \item{prob}{
Numeric scalar/vector; quantile probability between 0-1. The vector mode is only used in the specific case of \option{imager_func} = 'quantile' / 'quan' (with the speech marks). The scalar mode is used when \option{imager_func} = parquan. The latter only accepts a single value of \option{prob} but it is substantially faster to compute (using the \code{imager} package) and in general should be preferred.
}
  \item{dim_out}{
Integer vector; this defines the desired dimensions of the output image. If this is not provided then the output image is made to be the same size as the NAXIS1 and NAXIS2 arguments taken from \option{header_out} (which is usually what you will want TBH).
}
  \item{cores}{
Integer scalar; the number of cores to run top level stacking on. In general larger will be faster, assuming you have enough memory to hold that many images in memory. Once you are using around half of the available RAM you will see bigger performance loss due to memory pressure than gains due to more CPU power. Make sure not to set \option{cores} to more than the available cores, see \code{detectCores} in the \code{parallel} package (and probably no more than available cores minus 2). This needs to be used carefully in combination with the \option{cores_warp} argument since in total you do not want to set \option{cores} x \option{cores_warp} than you have have machine cores/threads available (probably machine cores/threads - 2). Default is half the available cores.
}
  \item{cores_warp}{
Integer scalar; the number of cores used for warping. This needs to be used carefully in combination with the \option{cores} argument since in total you do not want to set \option{cores} x \option{cores_warp} than you have have machine cores/threads available (probably machine cores/threads - 2).
}
  \item{multitype}{
Character scalar; options are 'fork' (often best for Unix like and non-Gui sessions) and 'cluster' (Windows, and when using a Gui like R-Studio or threaded task). 'cluster' is also meant to be a bit more stable (if you have issues) but there is more memory copying overhead, so will generally be slower.
}
  \item{Nbatch}{
Integer scalar; the size of processing batches to be used internally. This can be larger than \option{cores}, but usually it makes more sense to set \option{Nbatch} = \option{cores} since that will minimise memory and maximise CPU usage. If \option{return_all} = TRUE or \option{doclip} = TRUE then \option{Nbatch} will be set to the number of images (since all pre-stack projected images will either be returned or need to be processed for clipped stacking). The fastest possible speed (but maximal memory impact) will be \option{Nbatch} = \option{cores} = length(\option{image_list}). The slowest speed (but minimal memory impact) will be \option{Nbatch} = 1. A good compromise is to work out how many images will fit within half of the available RAM (Nsafe), and set \option{Nbatch} = \option{cores} = Nsafe (making sure not to set \option{cores} to more than the available cores, see \code{detectCores} in the \code{parallel} package).
}
  \item{keepcrop}{
Logical; this flag specifies if just the minimally cut down projection should be used internally. The only reason to set this to FALSE is if one of \option{return_all} or \option{dump_frames} is TRUE and the user once to inspect the warped frames easily (since they will be pixel aligned in this case). Setting \option{keepcrop} = TRUE will use much less memory and in general create stacks faster.
}
  \item{keep_extreme_pix}{
Logical; should a \option{cold} and \option{hot} pixel matrix be stored and returned? This can be useful for diagnosing and later masking bad cold and hot pixels from the stacked \option{image}. Potentially you can use these outputs to detect hot/cold biased pixels in the output \option{image} and replace them with the \option{cold}/\option{hot} pixel in the location, or alternatively just mask them entirely. Note when using \option{doclip} = TRUE then the \option{cold}/\option{hot} will be the coldest and hottest pixels remaining after the clipping is applied (so the originally most extreme pixels will no longer be present).
}
  \item{doclip}{
Logical; should the stack be cold/hot pixel clipped and restacked? Note \option{inVar_list} is required in this case, since the image variances are used to determine the clipping (see \option{clip_tol} for the selection applied). It only makes sense to attempt clipping when 3+ frames are being used, and in general more frames will provide better clipping.
}
  \item{clip_tol}{
Numeric scalar or vector; determines how should cold (\option{clip_tol}[1]) and hot (\option{clip_tol}[2]) pixels are clipped. The internal logic is (stack - cold).sqrt(inVar) > \option{clip_tol}[1] cold pixels and (hot - stack).sqrt(inVar) > \option{clip_tol}[2] hot pixels are clipped out for a subsequent re-stacking of the output \option{image}. If length is 1 then value is repeated for cold and hot pixels. In general you will want to set \option{clip_tol} in the range 5-100 since it is roughly speaking a sigma clip.
}
  \item{clip_dilate}{
Integer scalar; how much to dilate around cold/hot pixels you wish to mask out. This is only relevant when \option{doclip} is TRUE. The value of \option{clip_dilate} is the diameter of the circular kernal used, so it should be odd valued and >= 3 (which would be one pixel of dilation) to have any effect.
}
  \item{clip_sigma}{
Numeric scalar; clipping significance threshold for hot pixels. Only pixels where the \option{cold} pixel is below \option{clip_sigma} of the sky RMS will be flagged for hot pixel clipping. Generally you want this to be a value between 5-10, and it stops the cores of bright stars/galaxies being erroneously flagged as containing hot pixels. The trade off from setting this to higher values is you might start to lose the cores of bright objects, and the trade off setting this lower is you might not remove genuine hot pixels. In general this parameter should not need much fine-tuning.
}
  \item{return_all}{
Logical; should the pre stack projected images and inverse-variances be returns also? Note if this is TRUE then the stacking can no longer be computed in batch mode, so be sure that you have enough memory to hold all pre stack projected images and inverse-variances in memory (where note these will each have the dimensions of the desired output stacked image, which can be very large).
}
  \item{dump_frames}{
Logical; if TRUE then warped frames will be dumped out to the specified \option{dump_dir} as the stacking continues. This can be useful to inspect the pre-stack inputs (if too large for the \option{return_all} = TRUE option), or if the stack is failing unexpectedly. Frames will inherit the names image_warp_i, inVar_warp_i, weight_warp_i and exp_warp_i as appropriate (i.e. only if warping is occurring for these data). If \option{doclip} = TRUE then the boolean frame masks will also be dumped to mask_warp_i.
}
  \item{dump_dir}{
Character scalar; target directory to write out the dumped frames or (in the case of \code{propaneStackWarpMed}) read them in.
}
  \item{filelist}{
Character vector; vector of full paths of FITS files to analyse. This is written as an attribute (called \option{filename}) to the output 'Rfits_list' object.
}
  \item{dirlist}{
Character vector; vector of full paths of directories to analyse.
}
  \item{extlist}{
Integer vector; the extensions to use. If length 1 then it will be used for all \option{filelist}, but otherwise it should be the same length as the final \option{filelist}
}
  \item{pattern}{
Character vector; regular expressions to filter \option{filelist} by. Most people find it easier to work via \code{\link{glob2rx}} since then you can specify file wild cards in the usual shell way, e.g. glob2rx("*F200W*fits") becomes "^.*F200W.*fits$". Note an extra check that all files are FITS files (a pattern of ".fits$") will be made whether requested or not. Note you should nearly always use a leading and trailing '*' search because the pattern match is made on the full file path string, i.e. glob2rx('u*.fits') = '^u.*\\\\.fits$' would fail on '/path/to/file/u_GAMA.fits' (because of the '/path/to/file/' before 'u_GAMA.fits'), but glob2rx('*u*.fits') = '^.*u.*\\\\.fits$' would work as expected. For images dumped out by \code{\link{propaneStackWarpInVar}} when \option{dump_frames} = TRUE you will want to use glob2rx('*image_warp*') or glob2rx('*image_warp_*.fits').
}
  \item{recursive}{
Logical; if using \option{dirlist} should all sub-directories be checked recursively?
}
  \item{zap}{
Character vector; optional unique strings to zap out of the header. See \code{\link{Rfits_point}}.
}
  \item{chunk}{
Integer scalar; chunk sizes to maximally break the image up into to compute the median stacks. The sweet spot seems to be around 1000 - 2000 depending on sytem RAM and number of \option{cores}.
}
  \item{doweight}{
Logical; should weights (so number of frames contributing to pixel medians) be computed? If this is FALSE then the weight attribute of the list returns from \code{propaneStackWarpMed} will be NULL.
}
  \item{useCUTLO}{
Logical; should the \option{XCUTLO} and \option{YCUTLO} keyvalues in the header be used? This is usually what you want if using aligned WCS which have different cropping (like the \option{dump_frames} output of \code{\link{propaneStackWarpInVar}}), but if you have manually aligned the frames (e.g. for stacking stars to create a PSF) then you might want to ignore the \option{XCUTLO} and \option{YCUTLO} keyvalues created by \code{\link{Rfits}} methods.
}
  \item{\dots}{
Other arguments to pass into \code{\link{propaneWarp}}.
}
}
\details{
\code{propaneStackWarpInVar} is an inverse-variance weighted image stacking function, which has a lot of flexibility in terms of clipping and weighting.

The warping code used here is \code{\link{propaneWarp}}, and the stacking part is very similar to that used in \code{\link{propaneStackFlatInVar}}. It uses inverse-variance weighted stacking when \option{inVar_list} is provided, and a simple average stack otherwise.

\option{image_list} can be provided as a list of \code{Rfits_image} or \code{Rfits_pointer}. The advantage of the former is raw speed (no disk IO required) but the advantage of the latter is memory footprint (only the required pixels are loaded at any given stage of stacking).

\option{inVar_list}, \option{exp_list}, \option{weight_list} can all be provided as a scalar (repeated for all images) a vector of single values (used for all pixels of each image), a numeric/integer matrix (used for each image), an \code{Rfits_image} (the imDat part alone is used for each image) or an \code{Rfits_pointer} (the imDat part alone is used for each image).

\option{mask_list} can be provided as an integer matrix (used for each image), an \code{Rfits_image} (the imDat part alone is used for each image) or an \code{Rfits_pointer} (the imDat part alone is used for each image).

Which route is better to use depends on the memory available- passing in in-memory matrices will be faster than using an \code{Rfits_pointer}, but using a pointer means only the required pixels are loaded at a given time (reducing the memory footprint).

Care needs to be taken to choose appropriate values for \option{cores} and \option{Nbatch}. Ultimately at least a few copies of the full target image will need to be kept in memory at some point in this stacking process, since the stacking happens in RAM not on disk. For that reason there is an upper limit on how big the target image can be without hitting some pretty brutal performance issues (when you start using disk swap, things become very very slow). As an idea, a 10k x 10k pixel image with bitpix = -32 (sinlge precision) will take up 400 MB on disk and loaded in \code{R} 800 MB (since \code{R} only supports double precision floats internally). This means a 40k x 40k will require 12.8 GB of RAM. Given that number at least doubles during any \code{\link{propaneWarp}} process, and we possibly have a weight map and exposure to compute, you want to be pretty certain that x4 this value fits well within half of available RAM even if we are using \code{Nbatch} = 1. Even on a machine with 128 GB of RAM (my local machine) that is about as big an image as we could possibly target for stacking.

\code{propaneWarpDump} just caries out the frame warping and dumping part of \code{propaneStackWarpInVar}.

\code{propaneStackWarpMed} is a dedicated median combine stacking function. The difference in operation is that it only works on pre-dumped frames that will have been projected with \code{propaneStackWarpInVar} with the \option{dump_frames} option set to TRUE. The reason for this is that the median combine can generally only be run on quite small chunked subsets (usually in the range 1000 - 2000 and set with the \option{chunk} argument). Warping on the fly to such subsets is expensive (for reasonable sized inputs most pixels get thrown away each time), so it is faster to compute frame overlap on the dumped frames and combine them together as needed. The actual median combine is done using the \code{parmed} function of the \code{imager} package. If this has been compiled with OpenMP flags then you will also see massive speed up during the combining phase. The keen user will probably need to do some Googling to get this right, but for large stacking it could be worth the effort (for MacOS check out https://mac.r-project.org/openmp/).

\code{propaneStackWarpMed} can also be 'tricked' into stacking non-WCS aligned images. This might be useful when you are e.g. stacking stars that are on inherently different WCS, but you have already cutout and centred then. For this to work no \option{keyvalues_out} should be provided, and \option{useCUTLO} should be set to FALSE (just incase you have XCUTLO and YCUTLO keyvalues produced by \code{Rfits} in your header). The images to be stacked also need to be the same exact dimensions (as well as being centred sensibly on your source of interest, of course).

\code{propaneStackWarpFunc} is like a generic version of \code{propaneStackWarpMed}, where any of the ?imager.combine functions that collapse stacks of images down can be used rather than just \code{parmed} (which is the median function used in \code{propaneStackWarpMed}). For that reason, setting \option{imager_func} = imager::parmed will replicate exactly the behaviour of \code{propaneStackWarpMed}. Useful options to consider using are imager::parvar (variance) and imager::parsd (standard deviation). As a special case you can use \option{imager_func} = 'quantile' / 'quan' (with the speech marks) for quantiles (defined by \option{pro}). Also, \option{imager_func} = 'InVar' (with speech marks) will return the inverse-variance scaled correctly by exposures (which is what we use for general stacking). \option{imager_func} = 'waverage' will return a weighted average using the \option{weights} values for the per frame weights, but note this only works in regions where all frames contribute.
}
\value{
The output of \code{propaneStackWarpInVar} is an object of class 'ProPane' containing:

\item{image}{The stacked image, of class Rfits_image.}
\item{weight}{The number of available un-masked pixels in every output pixel of \option{image}. Class Rfits_image.}
\item{inVar}{If \option{inVar_list} is not NULL, the stacked inverse-variance of the output \option{image}. Class Rfits_image when present.}
\item{exp}{If \option{exp_list} is not NULL, the total exposure time of the output \option{image}. Class Rfits_image when present. Note this will be unmasked, so values of 0 means no photon ever hit the output \option{image}, i.e. masked pixels can have postive exposures if they sit within the detector domain.}
\item{cold}{If \option{keep_extreme_pix} = TRUE, the extreme hot pixel value for the input \option{image_list} that formed the final stacked \option{image}. If \option{keep_extreme_pix} = FALSE this will be NULL. Note when using \option{doclip} = TRUE then the \option{cold} will be the coldest pixels remaining after the clipping is applied (so the originally most extreme pixels will no longer be present). Class Rfits_image when present.}
\item{hot}{If \option{keep_extreme_pix} = TRUE, the extreme hot pixel value for the input \option{image_list} that formed the final stacked \option{image}. If \option{keep_extreme_pix} = FALSE this will be NULL. Note when using \option{doclip} = TRUE then the \option{hot} will be the hottest pixels remaining after the clipping is applied (so the originally most extreme pixels will no longer be present). Class Rfits_image when present.}
\item{clip}{The number of frames clipped from the stacked \option{image}. This is only relevant when \option{doclip} is TRUE. Values larger than 1 can indicate that both cold and hot pixels have been clipped in this location, and also that the dilation process has masked a pixel multiple times.}
\item{which_overlap}{Vector of IDs for the frames that overlap the target WCS.}
\item{time}{Numeric scalar; the total clock time taken (not collective CPU time) for the stacking.}
\item{Nim}{Integer scalar; number of images used in the stack.}
\item{dump_dir}{Character scalar; the input \option{dump_dir}.}

And optionally (if \option{returnall} = TRUE) these outputs are appended to the list:

\item{image_pre_stack}{The projected pre stack images.}
\item{inVar_pre_stack}{The projected pre stack inverse-variance.}
\item{exp_pre_stack}{The projected pre stack exposure times.}

The output of \code{propaneStackWarpMed} is an object of class 'ProPane' containing:

\item{image}{The median stacked image, of class Rfits_image.}
\item{weight}{The number of available un-masked pixels in every output pixel of \option{image}. Class Rfits_image.}
\item{which_overlap}{Vector of IDs for the frames that overlap the target WCS.}
\item{time}{Numeric scalar; the total clock time taken (not collective CPU time) for the stacking.}
\item{Nim}{Integer scalar; number of images used in the stack.}
\item{dump_dir}{Character scalar; the input \option{dump_dir}.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{propaneWarp}}, \code{\link{propaneStackFlat}}, \code{\link{propaneGenWCS}}
}
\examples{
\dontrun{
#Recommended to use Rfits (will be supported longer term):
library(Rfits)
library(Rwcs)

file_image = system.file('extdata', 'image.fits', package = "Rfits")
image = Rfits_read_image(file_image)

#Here we stack a cutout ontop of the original image. A silly use case, but shows how it
#works in basic terms.

stack = propaneStackWarpInVar(image_list=list(image, image[50:300,50:300]),
  magzero_out=0,
  direction='forward')

plot(stack$image) #this should look the same
plot(stack$weight, magmap=FALSE) #but centre of weight shows where WCS overlaps

#we can continue stacking the image, preserving the various objects correctly if
#they are passed in. Here we also cutout the original stack when re-stacking:

stack2 = propaneStackWarpInVar(image_list=list(stack$image, stack$image[20:350, 20:350]),
  weight_list=list(stack$weight$imDat,stack$weight$imDat[20:350, 20:350]),
  magzero_out=0,
  direction='forward')

plot(stack2$image)
plot(stack2$weight, magmap=FALSE)

#A better example with some real data that we can stack!

#Load in the target images

image_list = Rfits_make_list(dirlist = system.file('extdata/stack/', package="ProPane"),
                             extlist = 2) #extlist=2 because these are compressed images

#Define a distinct WCS to stack to (otherwise it uses the first in the image_list):

keyvalues_out = Rwcs_setkeyvalues(
    CRVAL1 = 36.8962,
    CRVAL2 = -5.1906,
    pixscale = 0.3,
    NAXIS1 = 2200,
    NAXIS2 = 2200
)

#Although see also propaneGenWCS for help doing this

stack = propaneStackWarpInVar(image_list = image_list,
                   inVar_list = 0.000468, #assume all have the same inVar
                   exp_list = 10, #all have the same exposure time of 10s
                   cores = 8,
                   keyvalues_out = keyvalues_out,
                   magzero_in = 30,
                   magzero_out = 23.9 #micro-jansky output
                   )

plot(stack$image, qdiff=TRUE)
plot(stack$inVar, magmap=FALSE)
plot(stack$weight, magmap=FALSE)
plot(stack$exp, magmap=FALSE)

#If we do not supply the keyvalues_out the target WCS will be the first image:

stack2 = propaneStackWarpInVar(image_list = image_list,
                   inVar_list = 0.000468, #assume all have the same inVar
                   exp_list = 10, #all have the same exposure time of 10s
                   cores = 8,
                   magzero_in = 30,
                   magzero_out = 23.9 #micro-jansky output
                   )

#Now we can easily see the depth improvment!

plot(image_list[[1]])
plot(stack2$image)

#We can compute the rough sky noise in the stacked images:

sd(stack$image$imDat[stack$image$imDat < 0.3], na.rm=TRUE) #stack 1
sd(stack2$image$imDat[stack2$image$imDat < 0.3], na.rm=TRUE) #stack 2

#Compare this to the stacked inverse-variance prediction:

1/sqrt(median(stack$inVar$imDat, na.rm=TRUE)) #stack 1
1/sqrt(median(stack2$inVar$imDat, na.rm=TRUE)) #stack 2

#Pretty similar! Note when changing pixel scales a lot this can be wrong due to pixel
#covariance, which we do not correct for.

#For median stacking we need to dump out the frames:

dump = propaneWarpDump(image_list = image_list,
                         cores = 8,
                         keyvalues_out = keyvalues_out,
                         magzero_in = 30,
                         magzero_out = 23.9 #micro-jansky output
)

stack_med = propaneStackWarpMed(dirlist=dump, keyvalues_out=keyvalues_out,
  pattern=glob2rx('*image_warp*'))

plot(stack_med$image, qdiff=TRUE)

#ProPane includes an example Swarp median stack with the same WCS projection:

swarp = Rfits_read_image(system.file('extdata/swarp_med.fits', package="ProPane"), ext=2)

plot(swarp, qdiff=TRUE)

# we can check the pixels are pretty similar!

magplot(stack_med$image$imDat[stack_med$image$imDat > 1],
  swarp$imDat[stack_med$image$imDat > 1], pch='.', log='xy', xlim=c(1,1e3),
  ylim=c(1,1e3), asp=1)
abline(0, 1, col='red')

maghist(stack_med$image$imDat[stack_med$image$imDat > 1] /
  swarp$imDat[stack_med$image$imDat > 1], xlim=5)

#Note Swarp tends to find slightly more flux (a couple percent) for the same pixels.
#This is because Swarp uses Lanczos interpolation by default, and imager bicubic.
#If this critically matters, consider upscaling your images with propaneRebin first!
#This works because simple upscaling uses a higher accuracy Lanczos interpolation.

#And we can measure the quantiles too:

stack_invar = propaneStackWarpFunc(dirlist=dump,
  keyvalues_out=keyvalues_out, imager_func='InVar', pattern=glob2rx('*image_warp*'))

#With a bit of thought, we can use this stack_invar in ProFound:
library(ProFound)

InVar = stack_invar$image$imDat
InVar = profoundImBlur(InVar) #Smooth out the InVar measurements a bit

#Run ProFound masking out regions where only 0/1 image contributes:

test_pro = profoundProFound(stack_mask, skyRMS=1/sqrt(InVar),
  mask=stack3$weight$imDat < 2, redosky=FALSE, magzero=23.9)

plot(test_pro)

#Now we are going to inject hot and cold pixels and try clipping:

image_list = Rfits_make_list(dirlist = system.file('extdata/stack/', package="ProPane"),
                              extlist = 2, pointer = FALSE) #load full images

set.seed(666) #set a random seed
selhot = sample(prod(dim(image_list[[1]])), 1e4)
image_list[[1]]$imDat[selhot] = 1e3 #inject hot pixels

selcold = sample(prod(dim(image_list[[2]])), 1e4)
image_list[[2]]$imDat[selcold] = -1e3 #inject cold pixels

plot(image_list[[1]][1e3,1e3,box=201], qdiff=TRUE)
plot(image_list[[2]][1e3,1e3,box=201], qdiff=TRUE)

#First we can see this will not look great without clipping:
stack3_noclip = propaneStackWarpInVar(image_list = image_list,
                   inVar_list = 0.000468, #assume all have the same inVar
                   exp_list = 10, #all have the same exposure time of 10s
                   cores = 8,
                   keyvalues_out = keyvalues_out,
                   magzero_in = 30,
                   magzero_out = 23.9, #micro-jansky output
                   keepcrop = FALSE, #so we can see the frames lined up
                   return_all = TRUE #so we keep the warped frames
                   )
plot(stack3_noclip$image[1e3,1e3,box=201], qdiff=TRUE)
plot(stack3_noclip$image_pre_stack[[1]][1e3,1e3,box=201], qdiff=TRUE) #hot frame
plot(stack3_noclip$image_pre_stack[[2]][1e3,1e3,box=201], qdiff=TRUE) #cold frame

#And much better with clipping:
stack3_clip = propaneStackWarpInVar(image_list = image_list,
                   inVar_list = 0.000468, #assume all have the same inVar
                   exp_list = 10, #all have the same exposure time of 10s
                   cores = 8,
                   keyvalues_out = keyvalues_out,
                   magzero_in = 30,
                   magzero_out = 23.9, #micro-jansky output
                   doclip = TRUE,
                   clip_tol=10
                   )

plot(stack3_clip$image[1e3,1e3,box=201], qdiff=TRUE)
plot(stack3_clip$clip[1e3,1e3,box=201])
}
}
\concept{ ~warp }
\concept{ ~stack }

