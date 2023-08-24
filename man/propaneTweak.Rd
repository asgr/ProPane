\name{propaneTweak}
\alias{propaneTweak}
\alias{propaneTweakImage}
\alias{propaneTweakCat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tweak the Translation Astrometry of an Image
}
\description{
Function to rapidly find translation and rotation tweaks directly in the image domain (i.e. it does not require a catalogue to be generated first). This is generally pretty fast, especially when running in integer shift mode.
}
\usage{
propaneTweak(image_ref, image_pre_fix, delta_max = c(3,0), quan_cut = 0.99, Nmeta = 3,
  WCS_match = TRUE, cores = 4, shift_int = TRUE, return_image = TRUE,
  direction = "backward", final_centre = TRUE, cutcheck = FALSE, quick = FALSE,
  qtype = 'optim', Niter = 1e4, verbose = TRUE)

propaneTweakImage(image_ref, image_pre_fix, delta_max = c(3,0), quan_cut = 0.99,
  Nmeta = 3, WCS_match = TRUE, cores = 4, shift_int = TRUE, return_image = TRUE,
  direction = "backward", final_centre = TRUE, cutcheck = FALSE, quick = FALSE,
  qtype = 'optim', Niter = 1e4, verbose = TRUE)

propaneTweakCat(cat_ref, cat_pre_fix, delta_max = c(100,0), mode = 'pix',
  keyvalues_ref = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_ref}{
Rfits_image / Numeric matrix; static reference image (must have exactly the same dimensions as \option{image_pre_fix} if \option{shift_int} or \option{WCS_match} is TRUE). If \option{WCS_match} is FALSE then \option{image_ref} must be a Rfits_image.
}
  \item{image_pre_fix}{
Rfits_image / Numeric matrix; image to be tweaked (must have exactly the same dimensions as \option{image_ref} if \option{shift_int} or \option{WCS_match} is TRUE). If \option{WCS_match} is FALSE then \option{image_pre_fix} must be a Rfits_image.
}
  \item{cat_ref}{
Numeric matrix / data.frame; simple two column reference catalogue containing x/RA [,1] and y/Dec [,2] columns. xy corresponds to \option{mode} = 'pix', 'RA/Dec to  \option{mode} = 'coord'. The \option{mode} type must be the same as \option{cat_pre_fix}. If \option{mode} = 'coord' then \option{keyvalues_ref} must also be provided.
}
  \item{cat_pre_fix}{
Numeric matrix / data.frame; simple two column tweak catalogue containing x/RA [,1] and y/Dec [,2] columns. xy corresponds to \option{mode} = 'pix', 'RA/Dec to  \option{mode} = 'coord'. The \option{mode} type must be the same as \option{cat_ref}. If \option{mode} = 'coord' then \option{keyvalues_ref} must also be provided.
}
  \item{delta_max}{
Numeric vector; the first element is maximum number of pixels to explore in +/-x and +/- y; the second element is the maximum rotation to apply in degrees. If \option{delta_max}[2] = 0 then only translation optimisation is computed (no rotation).
}
  \item{quan_cut}{
Numeric scalar/vector; when length 1 a bright pixel thesholding to apply to \option{image_ref} and \option{image_pre_fix}. Generally you want this quite high so that only cores of bright objects are compared. Suggested range to try (depending on image crowding) is 0.9 - 0.999 (see also \option{cutcheck} to aid this process). If length 2 then the first element (\option{quan_cut}[1]) is the bright pixel thresholding, and the second element (\option{quan_cut}[2]) is used to clip artificially bright pixels. Generally the seconed value should be in the range 0.9999 - 0.99999.
}
  \item{Nmeta}{
Integer scalar; how many shifts of the search box can be applied. This is only relevant when \option{shift_int} is TRUE.
}
  \item{WCS_match}{
Logical; specifies whether we are expecting the input WCS to be matched between \option{image_ref} and \option{image_pre_fix} or not. If \option{WCS_match} is TRUE then the image itself is directly manipulated and the WCS is left alone. If \option{WCS_match} is FALSE then the WCS header is modified and \option{image_pre_fix} is projected onto the WCS of \option{image_ref} for each iteration (therefore \option{image_ref} and \option{image_pre_fix} must both be Rfits_image objects in this case). \option{WCS_match} = FALSE is often the more useful way to tweak images (since often images are known not to have exactly matching WCS because they are dithered etc), but it is quite a bit slower (about x5 longer to run) since the projection has to go through \code{\link{propaneWarp}} and therefore \code{wcslib}.
}
  \item{cores}{
Integer scalar; the number of cores to run on.
}
  \item{shift_int}{
Logical; should the shift applied be purely integer (TRUE) or fractional pixel (FALSE). Integer shifts are much faster to compute, and are better at conserving flux and not introducing pixel covariance, so think carefully before setting to FALSE.
}
  \item{return_image}{
Logical; should the tweaked image be returned as an output? You might only want the solution (for applying to another process), in which case set to FALSE for speed. The returned image will either be pixel shifted if \option{image_pre_fix} is a numeric matrix, or will have its header modified if \option{image_pre_fix} is a Rfits_image using \code{\link{propaneWCSmod}}.
}
  \item{direction}{
Character scalar; "forward" or "backward", see \code{imwarp}.
}
  \item{final_centre}{
Logical; should a final round of costs be computed once the best tweak is found. This is only relevant if \option{shift_int} = TRUE, and can be useful if you want to explore the output \option{cost_mat}.
}
  \item{cutcheck}{
Logical; if TRUE then quick plots showing the results of \option{quancut} are displayed. If the good quality sources to WCS match against are not present then consider lowering \option{quancut} (suggestion of 0.95 of 0.9). If a lot of noise and artefacts are present then consider increasing it. Note if this is TRUE then the function will exit after displaying the plots returning NULL (i.e. no tweaking is actually done).
}
  \item{quick}{
Logical; determines whether an approximate 'quick' solution is found. This will often be very good, but it does not compute a full image warp (it rounds pixels to the nearest match).
}
  \item{qtype}{

}
  \item{Niter}{
Integer scalar; if \option{quick} = TRUE, then this sets the maximum number of iterations to try. If it converges faster then the \option{Niter} will not be reached, so usually best not to reduce it from the default.
}
  \item{mode}{
Logical; the type of input for \option{cat_ref} and \option{cat_pre_fix}, either 'pix' (for x/y pixels) or 'coord' (For RA/Dec in degrees). Both catalogues must have the same data \option{mode}. If \option{mode} = 'coord' then \option{keyvalues_ref} must also be provided.
}
 \item{keyvalues_ref}{
List; header values that detail the WCS for \option{cat_ref}. This must be provided when \option{mode} = 'coord'.
}
  \item{verbose}{
Logical; should verbose outputs be printed as it runs?
}
  \item{\dots}{
Other arguments to pass to \code{coordmatch}.
}
}
\details{
Image tweaking is useful before stacking with e.g. \code{\link{propaneStackWarpInVar}} if images are slightly misaligned in translation, but not rotated or distorted otherwise. This is quite common with HST and JWST images when guide starts are switched between observing blocks.

If using the output of this to fix the header CRPIX then the fixed/new CRPIX belonging to \option{image_pre_fix} becomes \option{image_pre_fix}$keyvalues$CRPIX[1-2] (fixed) = \option{image_pre_fix}$keyvalues$CRPIX[1-2] (old) - tweak$optim_out$par[1-2] (i.e. we subtract the tweaked solution). In this case you might need to replace the raw list item with either the output of Rfits_keyvalues_to_header(Rfits_header_to_raw(keyvalues_fixed)) or set to NULL (this is because, if present, raw header items take precedence).
}
\value{
For \code{propaneTweak} and \code{propaneTweakImage} (which are the same function):

When \option{shift_int} is TRUE:

\item{optim_out}{The output of the actual optimisation (includes the best par etc). This can be used to update the relevent image header of \option{image_pre_fix} via \code{\link{propaneWCSmod}} (see Examples)}
\item{image_post_fix}{The tweaked image (if \option{return_image} is TRUE). This will either have the image pixels directly manipulated or the WCS header, depending on whether a Rfits_image was passed in for \option{image_pre_fix}. There are 3 potential outputs: if the input \option{image_pre_fix} was an image matrix then obviously the output has to be a modified image (there is no header), but if the input \option{image_pre_fix} is an Rfits_image then \option{WCS_match} = FALSE means the header will be changed, and \option{WCS_match} = TRUE means the image pixels will be adjusted direclty.}
\item{time}{The total running time in seconds.}
\item{cost_mat}{The cost outputs for the final iteration. Thie has dimensions x = delta_x, y = delta_y, z = cost.}
\item{at_lim}{Logical; has the solution hit a limit. If so maybe re-run with larger \option{delta_max} and/or larger \option{Nmeta}.}

When \option{shift_int} is FALSE:

\item{optim_out}{The output of the actual optimisation (includes the best par etc).}
\item{image_post_fix}{The tweaked image (if \option{return_image} is TRUE).}
\item{time}{The total running time in seconds.}

For \code{propaneTweakCat}:

The output of the actual optimisation (includes the best par). This can be used to update the relevent image header of \option{cat_pre_fix} via \code{\link{propaneWCSmod}} (see Examples).
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{propaneTran}}
}
\examples{
library(Rfits)
library(magicaxis)

image = Rfits_read_image(system.file('extdata/VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))
plot(image)

# below we shift the image by +2, +2 and rotate it by 5 degrees:
image_tran = image
image_tran$imDat = propaneTran(image_tran$imDat, 2, 2, 5)
plot(image_tran)

\dontrun{
# Now we tweak the solution (this will take about 20 seconds):

tweak = propaneTweak(image, image_tran, delta_max=c(5,10), shift_int=FALSE)

# Our tweaked image is pretty close to this:

best = tweak$optim_out$par
print(best)

# We can use the output to manipulate the pixels:

image_fix = image
image_fix$imDat = propaneTran(image_tran$imDat, best[1], best[2], best[3])
plot(image_fix)

# We can also modify the WCS rather than the image pixels:

image_fix2 = propaneWCSmod(image_tran, best[1], best[2], best[3]) #correct WCS
plot(image_fix2)

# In the last example note how the WCS grid follows the original sources correctly

# To make all of the above easier we can request propaneTweak to directly modify our
# WCS passed in (which is the default):

plot(tweak$image_post_fix)

# Now we will look at a case where the pixels are matched but the WCS is not:

image_tran2 = propaneWCSmod(image, 2, 2, 5)
plot(image_tran2)

# Note we set WCS_match to FALSE since it is the WCS itself that is wrong.
# The tweak of the WCS directly will takes longer than before (about 30 seconds):

tweak2 = propaneTweak(image, image_tran2, delta_max=c(5,10), shift_int=FALSE,
  WCS_match=FALSE)

# The answer almost perfectly corrects are adjustment to the WCS:

print(tweak2$optim_out$par)

# And the image and WCS looks very nice now:

plot(tweak2$image_post_fix)

# we can compare this to the 'quick' tweaking method:

tweak3 = propaneTweak(image, image_tran2, delta_max=c(5,10), shift_int=FALSE,
  WCS_match=FALSE, quick=TRUE)

print(tweak3$optim_out$par)

plot(tweak3$image_post_fix)

# Now some catalogue level tweaking:

library(ProFound)

plot(image_tran)

pro_im = profoundProFound(image)$segstats
pro_im_tran = profoundProFound(image_tran)$segstats

tweak4 = propaneTweakCat(pro_im[,c("xcen", "ycen")], pro_im_tran[,c("xcen", "ycen")],
  delta_max = c(5,20))
best3 = tweak4$par
best3
# The above is very fast, but obviously a bit less accurate (there is less information
# and signal in a catalogue versus the full image).

image_fix3 = image
image_fix3$imDat = propaneTran(image_tran$imDat, best3[1], best3[2], best3[3])
plot(image_fix)

# Or also:

image_fix4 = propaneWCSmod(image_tran, best3[1], best3[2], best3[3]) #correct WCS
plot(image_fix4)

# Where the WCS is now rotated to follow that of the original image.
}
}

\concept{ ~tweak }

