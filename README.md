# ProPane

## Synopsis

Package of stacking related functions, replacing older **Rwcs** and **ProFound** related functions.

The old to new function conversions are:

*Rwcs_stack* &rarr; *propaneStackWarpInVar*

*Rwcs_stack_median* &rarr; *propaneStackWarpMed*

*profoundMakeStack* &rarr; *propaneStackFlatInVar*

*profoundCombine* &rarr; *propaneStackFlatFunc*

AS of 06/02/2023 the functions above are identical to **Rwcs** v1.6.0 and **ProFound** v1.20.4, but from that point only the **ProPane** versions of these functions will be maintained. Much longer term (probably the next moderate version increments) these older function will be removed from **Rwcs** and **ProFound** respectively.

## Installation

You will need a number of other packages for the warping and stacking functions to work, in particular Rwcs and Rfits (on asgr GitHub) and imager (on CRAN but also on asgr GitHub).
