# ProPane

**ProPane** is a general purpose image warping and stacking package. It uses the **wcslib** **C** library for projections (all legal ones are supported) via the **Rwcs** package, and uses the threaded **Cimg** **C++** library via the **imager** library to do image warping. Most of the CPU expensive routines are parallelised at the **R** level, and some aspects are threaded at the **C** level via **OpenMP**. **imager** in particular will need to be built with **OpenMP** threading support to make use of these speed-ups, e.g. for MacOS you will want to look at this page https://mac.r-project.org/openmp/.

## Synopsis

Package of warping and stacking related functions. The initial version of this package involved converting older **Rwcs** and **ProFound** related functions. The non-**ProPane** versions of these functions will eventually be deprecated and removed, and from 06/02/2023 only the newer **ProPane** versions will be maintained (bug fixed etc) and enhanced.

The old to new function conversions are:

*Rwcs_warp* → *propaneWarp*

*Rwcs_stack* → *propaneStackWarpInVar*

*Rwcs_rebin* → *propaneRebin*

*Rwcs_tran* → *propaneTran*

*Rwcs_tweak* → *propaneTweak*

*Rwcs_warp_promo* → *propaneWarpProPane*

*Rwcs_stack_median* → *propaneStackWarpMed*

*profoundMakeStack* → *propaneStackFlatInVar*

*profoundCombine* → *propaneStackFlatFunc*

As of 06/02/2023 the functions above are identical to **Rwcs** v1.6.0 and **ProFound** v1.20.4, but from that point only the **ProPane** versions of these functions will be maintained. Much longer term (probably the next moderate version increments) these older function will be removed from **Rwcs** and **ProFound** respectively. Notably all the **Rwcs** equivalent functions refer to 'ProMo' class stacked objects, whereas in **ProPane** the new objects are 'ProPane' class stacked objects.

## Installation

You will need a number of other packages for the warping and stacking functions to work, in particular **Rwcs** and **Rfits** (on asgr GitHub) and **imager** (on CRAN but also on asgr GitHub). Once you have these you will be able to install **ProPane** with the following:

``` r
# install.packages("remotes")
remotes::install_github("asgr/ProPane")
```
## Example

Load in the target images:

```r
image_list = Rfits_make_list(dirlist = system.file('extdata/stack/', package="ProPane"),
                             extlist = 2) #extlist=2 because these are compressed images
```

Define a distinct WCS to stack to (otherwise it uses the first in the image_list):

```r
keyvalues_out = Rwcs_setkeyvalues(
    CRVAL1 = 36.8962,
    CRVAL2 = -5.1906,
    pixscale = 0.3,
    NAXIS1 = 2200,
    NAXIS2 = 2200
)
```

Run a parallel inverse variance weighted stack (this will take a few seconds):

```r
stack = propaneStackWarpInVar(image_list = image_list,
                   inVar_list = 0.000468, #assume all have the same inVar
                   exp_list = 10, #all have the same exposure time of 10s
                   cores = 8,
                   keyvalues_out = keyvalues_out,
                   magzero_in = 30,
                   magzero_out = 23.9 #micro-jansky output
                   )
```

Make some useful plots of the output:

```r
plot(stack$image, qdiff=TRUE)
plot(stack$inVar, magmap=FALSE)
plot(stack$weight, magmap=FALSE)
plot(stack$exp, magmap=FALSE)
```

For more complicated use cases please see the extended vignette at https://rpubs.com/asgr/999429.
