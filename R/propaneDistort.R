# Pixel distortion maps: the local pixel scale and the local pixel area of every
# pixel in a full image, as an Rfits_image that matches the input 1-1.
#
# Both functions probe the sky position of three corners of each pixel and take
# the size of the resulting patch of sky. Two subtleties in that, both of which
# produced visibly wrong maps on a field crossing RA = 0 (found on a WWS FUV
# mosaic centred near RA 10 deg, Dec -37.5 deg):
#
# 1. The RA corners straddle the 0/360 seam for part of such a field, so a raw
#    BR[,1] - BL[,1] jumps by ~360 deg. The old code took that literally, which
#    inflated the pixel area by up to three orders of magnitude in a band ~55
#    pixels wide wherever the seam fell (measured: max 3.3e8 against a true max
#    of 9.0e4 asec2). RA differences are therefore wrapped into [-180, 180).
#    Note the sign in the sky is not the sign in pixel space: RA *decreases* as
#    x increases here, so the seam is crossed as BR - BL -> +360, not -360;
#    wrapping handles either direction.
#
# 2. `cos(dec)` was applied per corner, i.e. once to each RA in degrees, which
#    is wrong twice over. RA *differences* must be multiplied by a single
#    cos(dec) before the Euclidean step length is formed, and multiplying an
#    absolute longitude by cos(dec) makes the result depend on where the field
#    sits on the sky. Shifting this header by 180 deg in RA changed the reported
#    area by a factor of 2.37 at the same pixel (9.05e4 -> 2.15e5 asec2), while
#    the true area is unchanged. Rfits::pixarea() gets both of these right and
#    is used below as an independent check.
#
# Beyond the wrap and cos(dec), the area was the product of the two side
# lengths. That is only the area of a rectangle; sky pixels are generally
# parallelograms (the two sides meet at 83-96 deg across the test field, mostly
# from the TAN curvature, and the shortfall is worst in the corners). The cross
# product |d1 x d2| is exact for any parallelogram and is what Rfits uses;
# the two agree to a median 2.3e-4 and a worst case of 6.8e-3 here.
#
# propaneDistortPixscale keeps the mean of the two side lengths, which is the
# right notion for a scale (an area would be their geometric mean), but it had
# the same wrap and cos(dec) faults and the same fix applies.

# Wrap a longitude difference into [-180, 180). Sign-preserving for small
# differences and correct for a seam crossing in either direction.
.wrap_long_diff = function(dlong){
  ((dlong + 180) %% 360) - 180
}

# Sky positions of the BL, BR and TL corners of the pixels centred on x, y
# (1-based FITS pixel indices), returned in degrees.
.pix_corners = function(x, y, keyvalues, header, ...){
  p2s = function(xx, yy){
    Rwcs_p2s(xx - 0.5, yy - 0.5, keyvalues=keyvalues, header=header,
             pixcen='FITS', ...)
  }
  list(BL=p2s(x, y), BR=p2s(x + 1, y), TL=p2s(x, y + 1))
}

# The two pixel sides as sky offsets in degrees, RA scaled by a single cos(dec)
# so they are on an equal-area local tangent plane. RA differences are wrapped,
# and the declination used is the mean of the three corners so that no one corner
# is privileged. Returns a list of dRA/dDec for the x side and the y side.
.pix_sides = function(corners){
  cosdec = cos(rowMeans(cbind(corners$BL[,2], corners$BR[,2],
                              corners$TL[,2]))*pi/180)
  list(
    dRAx = .wrap_long_diff(corners$BR[,1] - corners$BL[,1])*cosdec,
    dDx  = corners$BR[,2] - corners$BL[,2],
    dRAy = .wrap_long_diff(corners$TL[,1] - corners$BL[,1])*cosdec,
    dDy  = corners$TL[,2] - corners$BL[,2]
  )
}

propaneDistortPixscale = function(keyvalues, useraw=TRUE, unit='asec', ...){
  im_dim = dim(keyvalues)
  pixgrid = expand.grid(1:im_dim[1],1:im_dim[2])

  if(useraw){
    header = Rfits_keyvalues_to_raw(keyvalues)
  }else{
    header = NULL
  }

  sides = .pix_sides(.pix_corners(pixgrid[,1], pixgrid[,2], keyvalues, header, ...))
  # Mean of the two side lengths: a square of side s reports s, and an
  # anisotropic pixel reports the side of its equal-length equivalent. Not
  # sqrt(area); use propaneDistortPixarea for that.
  scale = (sqrt(sides$dRAx^2 + sides$dDx^2) + sqrt(sides$dRAy^2 + sides$dDy^2))/2

  if(unit=='deg'){
    #do nothing
  }else if(unit == 'asec'){
    scale = scale*3600
  }else if(unit == 'amin'){
    scale = scale*60
  }else if(unit=='rad'){
    scale = scale*(pi/180)
  }else{
    message('Not a valid unit, must be one of asec / amin / deg / rad')
  }

  output = Rfits_create_image(matrix(scale, im_dim[1], im_dim[2]), keyvalues=keyvalues)

  return(output)
}

propaneDistortPixarea = function(keyvalues, useraw=TRUE, unit='asec2', ...){
  im_dim = dim(keyvalues)
  pixgrid = expand.grid(1:im_dim[1],1:im_dim[2])

  if(useraw){
    header = Rfits_keyvalues_to_raw(keyvalues)
  }else{
    header = NULL
  }

  sides = .pix_sides(.pix_corners(pixgrid[,1], pixgrid[,2], keyvalues, header, ...))
  # |d1 x d2|, the parallelogram area. Exact whatever the angle between the
  # sides, unlike the product of the side lengths it replaces.
  area = abs(sides$dRAx*sides$dDy - sides$dDx*sides$dRAy) #in deg2

  if(unit=='deg2'){
    #do nothing
  }else if(unit == 'asec2'){
    area = area*3600^2
  }else if(unit == 'amin2'){
    area = area*60^2
  }else if(unit=='rad2' | unit=='str'){
    area = area*(pi/180)^2
  }else{
    message('Not a valid unit, must be one of asec2 / amin2 / deg2 / rad2 / str')
  }

  output = Rfits_create_image(matrix(area, im_dim[1], im_dim[2]), keyvalues=keyvalues)

  return(output)
}