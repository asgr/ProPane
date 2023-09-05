propaneCropOverlap = function(image_ref, image_warp, ...){
  image_ref = Rfits_crop(image_ref)
  image_warp = propaneWarp(image_warp, image_ref$keyvalues, keepcrop=TRUE, extratight=TRUE, ...)

  xlo = image_warp$keyvalues$XCUTLO
  xhi = image_warp$keyvalues$XCUTHI
  ylo = image_warp$keyvalues$YCUTLO
  yhi = image_warp$keyvalues$YCUTHI

  if(xlo != 1L | xhi != dim(image_ref)[1] | ylo != 1L & yhi != dim(image_ref)[2]){
    image_ref = image_ref[xlo:xhi, ylo:yhi]
  }

  return(list(image_ref=image_ref, image_warp=image_warp))
}
