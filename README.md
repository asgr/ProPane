# ProPane

## Synopsis

Package of warping and stacking related functions, converting (and eventually replacing) older **Rwcs** and **ProFound** related functions.

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

As of 07/02/2023 the functions above are identical to **Rwcs** v1.6.0 and **ProFound** v1.20.4, but from that point only the **ProPane** versions of these functions will be maintained. Much longer term (probably the next moderate version increments) these older function will be removed from **Rwcs** and **ProFound** respectively. Notably all the **Rwcs** equivalent functions refer to 'ProMo' class stacked objects, whereas in **ProPane** the new objects are 'ProPane' class stacked objects.

## Installation

You will need a number of other packages for the warping and stacking functions to work, in particular Rwcs and Rfits (on asgr GitHub) and imager (on CRAN but also on asgr GitHub).
