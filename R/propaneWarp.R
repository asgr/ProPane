.warpfunc_in2out = function(x, y, header_in=NULL, WCSref_in=NULL, header_out=NULL, WCSref_out=NULL, cores=1) {
  radectemp = Rwcs_p2s(x, y, header = header_in, WCSref = WCSref_in, cores = cores)
  xy_out = Rwcs_s2p(radectemp[,1], radectemp[,2], header = header_out, WCSref = WCSref_out, cores = cores)
  return(xy_out)
}
.warpfunc_out2in = function(x, y, header_in=NULL, WCSref_in=NULL, header_out=NULL, WCSref_out=NULL, cores=1) {
  radectemp = Rwcs_p2s(x, y, header = header_out, WCSref = WCSref_in, cores = cores)
  xy_out = Rwcs_s2p(radectemp[,1], radectemp[,2], header = header_in, WCSref = WCSref_out, cores = cores)
  return(xy_out)
}

propaneWarp = function(image_in, keyvalues_out=NULL, header_out = NULL, dim_out = NULL,
                       direction = "auto", boundary = "dirichlet", interpolation = "cubic",
                       doscale = TRUE, dofinenorm = TRUE, plot = FALSE, dotightcrop = TRUE,
                       keepcrop = FALSE, extratight = FALSE, WCSref_out = NULL, WCSref_in = NULL,
                       magzero_out = NULL, magzero_in = NULL, blank=NA, warpfield=NULL,
                       warpfield_return=FALSE, cores=1, checkWCSequal=FALSE, ...)
{
  if(!requireNamespace("Rwcs", quietly = TRUE)){
    stop("The Rwcs package is needed for this function to work. Please install it from GitHub asgr/Rwcs", call. = FALSE)
  }

  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("The imager package is needed for this function to work. Please install it from CRAN.", call. = FALSE)
  }

  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop("The Rfits package is needed for this function to work. Please install it from GitHub asgr/Rfits", call. = FALSE)
  }

  if(! inherits(image_in, c('Rfits_image', 'Rfits_pointer'))){
    stop('image_in must be either Rfits_image or Rfits_pointer!')
  }

  keyvalues_in = image_in$keyvalues
  keyvalues_in = keyvalues_in[!is.na(keyvalues_in)]
  header_in = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_in))

  if(is.character(header_out)){
    if(!is.null(keyvalues_out)){
      message('Using header_out and ignoring keyvalues_out!')
    }
    if(length(header_out) == 1){
      keyvalues_out = Rfits_header_to_keyvalues(Rfits_raw_to_header(header_out))
    }else if(length(header_out) > 1){
      keyvalues_out = Rfits_hdr_to_keyvalues(header_out)
    }
  }

  if(is.null(keyvalues_out) & is.null(header_out)){
    keyvalues_out = options()$current_keyvalues
    header_out = options()$current_header
  }

  keyvalues_out = keyvalues_out[!is.na(keyvalues_out)]
  header_out = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_out))

  if(checkWCSequal){
    if(Rfits_key_match(keyvalues_out, keyvalues_in,
                       check = c('NAXIS1',
                                 'NAXIS2',
                                 'CRPIX1',
                                 'CRPIX2',
                                 'CRVAL1',
                                 'CRVAL2',
                                 'CTYPE1',
                                 'CTYPE2',
                                 'CUNIT1',
                                 'CUNIT1',
                                 'CD1_1',
                                 'CD1_2',
                                 'CD2_1',
                                 'CD2_2'
                        )
                      )
    ){
      message('WCS appears to be the same, directly returning input!')
      image_in$keyvalues$XCUTLO = 1L
      image_in$keyvalues$XCUTHI = dim(image_in)[1]
      image_in$keyvalues$YCUTLO = 1L
      image_in$keyvalues$YCUTHI = dim(image_in)[2]

      image_in$hdr = Rfits_keyvalues_to_hdr(image_in$keyvalues)
      image_in$header = Rfits_keyvalues_to_header(image_in$keyvalues)
      image_in$raw = Rfits_header_to_raw(Rfits_keyvalues_to_header(image_in$keyvalues))

      image_in$keynames = names(image_in$keyvalues)

      image_in$keycomments$XCUTLO = 'Low image x range'
      image_in$keycomments$XCUTHI = 'High image x range'
      image_in$keycomments$YCUTLO = 'Low image y range'
      image_in$keycomments$YCUTHI = 'High image y range'

      return(invisible(image_in))
    }
  }

  if (!is.null(keyvalues_out) & is.null(dim_out)){
    if(is.null(keyvalues_out$ZNAXIS1)){
      NAXIS1 = keyvalues_out$NAXIS1
    }else{
      NAXIS1 = keyvalues_out$ZNAXIS1
    }
    if(is.null(keyvalues_out$ZNAXIS2)){
      NAXIS2 = keyvalues_out$NAXIS2
    }else{
      NAXIS2 = keyvalues_out$ZNAXIS2
    }
    dim_out = c(NAXIS1, NAXIS2)
  }else{
    keyvalues_out$NAXIS1 = dim_out[1]
    keyvalues_out$NAXIS2 = dim_out[2]
  }

  if(is.null(dim_out)){
    stop('Missing NAXIS1 / NAXIS2 in header keyvalues! Specify dim_out.')
  }

  if(dotightcrop){
    suppressMessages({
      BL_out = Rwcs_p2s(0, 0, header=header_out, pixcen='R', WCSref=WCSref_out)
      TL_out = Rwcs_p2s(0, dim_out[2], header=header_out, pixcen='R', WCSref=WCSref_out)
      TR_out = Rwcs_p2s(dim_out[1], dim_out[2], header=header_out, pixcen='R', WCSref=WCSref_out)
      BR_out = Rwcs_p2s(dim_out[1], 0, header=header_out, pixcen='R', WCSref=WCSref_out)
    })

    corners_out = rbind(BL_out, TL_out, TR_out, BR_out)
    tightcrop_out = ceiling(Rwcs_s2p(corners_out, header=header_in, pixcen='R', WCSref=WCSref_in))

    min_x_out = max(1L, min(tightcrop_out[,1]))
    max_x_out = min(dim(image_in)[1], max(tightcrop_out[,1]))
    min_y_out = max(1L, min(tightcrop_out[,2]))
    max_y_out = min(dim(image_in)[2], max(tightcrop_out[,2]))

    if(min_x_out != 1 | max_x_out != dim(image_in)[1] | min_y_out != 1 | max_y_out != dim(image_in)[2]){
      if(inherits(image_in, 'Rfits_pointer')){
        image_in = image_in[c(min_x_out, max_x_out), c(min_y_out, max_y_out), header=TRUE]
      }else{
        image_in = image_in[c(min_x_out, max_x_out), c(min_y_out, max_y_out)]
      }

      keyvalues_in = image_in$keyvalues
      header_in = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_in))
    }else{
      if(inherits(image_in, 'Rfits_pointer')){
        image_in = image_in[,]
      }
    }

    suppressMessages({
      BL_in = Rwcs_p2s(0, 0,header=header_in, pixcen='R', WCSref=WCSref_in)
      TL_in = Rwcs_p2s(0, dim(image_in)[2], header=header_in, pixcen='R', WCSref=WCSref_in)
      TR_in = Rwcs_p2s(dim(image_in)[1], dim(image_in)[2], header=header_in, pixcen='R', WCSref=WCSref_in)
      BR_in = Rwcs_p2s(dim(image_in)[1], 0, header=header_in, pixcen='R', WCSref=WCSref_in)
    })

    corners_in = rbind(BL_in, TL_in, TR_in, BR_in)
    tightcrop_in = ceiling(Rwcs_s2p(corners_in, header=header_out, pixcen='R', WCSref=WCSref_out))

    min_x_in = max(1L, min(tightcrop_in[,1]))
    max_x_in = max(min_x_in + dim(image_in)[1] - 1L, range(tightcrop_in[,1])[2])
    min_y_in = max(1L, min(tightcrop_in[,2]))
    max_y_in = max(min_y_in + dim(image_in)[2] - 1L, range(tightcrop_in[,2])[2])

    # new code should be more efficient!

    if(!isTRUE(keyvalues_out$ZIMAGE)){
      keyvalues_out$NAXIS1 = max_x_in - min_x_in + 1L
      keyvalues_out$NAXIS2 = max_y_in - min_y_in + 1L
    }else{
      keyvalues_out$ZNAXIS1 = max_x_in - min_x_in + 1L
      keyvalues_out$ZNAXIS2 = max_y_in - min_y_in + 1L
    }

    keyvalues_out$CRPIX1 = keyvalues_out$CRPIX1 - min_x_in + 1L
    keyvalues_out$CRPIX2 = keyvalues_out$CRPIX2 - min_y_in + 1L

    header_out = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_out))

    image_out = list(
      imDat = matrix(c(blank,image_in$imDat[0]), max_x_in - min_x_in + 1L, max_y_in - min_y_in + 1L),
      keyvalues = keyvalues_out,
      hdr = Rfits_keyvalues_to_hdr(keyvalues_out),
      header = Rfits_keyvalues_to_header(keyvalues_out),
      raw = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_out)),
      keynames = names(keyvalues_out),
      keycomments = as.list(rep('', length(keyvalues_out)))
    )
    names(image_out$keycomments) = image_out$keynames
    class(image_out) = c('Rfits_image', class(image_out))
  }else{
    if(inherits(image_in, 'Rfits_pointer')){
      image_in = image_in[,]
    }

    image_out = list(
      imDat = matrix(c(blank,image_in$imDat[0]), max(dim(image_in)[1], dim_out[1]), max(dim(image_in)[2], dim_out[2])),
      keyvalues = keyvalues_out,
      hdr = Rfits_keyvalues_to_hdr(keyvalues_out),
      header = Rfits_keyvalues_to_header(keyvalues_out),
      raw = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_out)),
      keynames = names(keyvalues_out),
      keycomments = as.list(rep('', length(keyvalues_out)))
    )
    names(image_out$keycomments) = image_out$keynames
    class(image_out) = c('Rfits_image', class(image_out))

    min_x_in = 1L
    max_x_in = dim_out[1]
    min_y_in = 1L
    max_y_in = dim_out[2]
  }

  dim_min_x_in = min(dim(image_in)[1], dim(image_out$imDat)[1])
  dim_min_y_in = min(dim(image_in)[2], dim(image_out$imDat)[2])

  image_out$imDat[1:dim_min_x_in, 1:dim_min_y_in] = image_in$imDat[1:dim_min_x_in, 1:dim_min_y_in]
  rm(image_in)

  if(!is.null(magzero_in) & !is.null(magzero_out)){
    image_out$imDat = image_out$imDat*10^(-0.4*(magzero_in - magzero_out))
    keyvalues_out$MAGZERO = magzero_out
  }

  suppressMessages({
    pixscale_in = Rwcs_pixscale(keyvalues=keyvalues_in)
    pixscale_out = Rwcs_pixscale(keyvalues=keyvalues_out)
  })

  if (direction == "auto") {
    if (pixscale_in < pixscale_out) {
      direction = "forward"
    }
    if (pixscale_in >= pixscale_out) {
      direction = "backward"
    }
  }

  if(is.null(warpfield)){
    pix_grid = expand.grid(1:dim(image_out$imDat)[1], 1:dim(image_out$imDat)[2])

    if (direction == "forward") {
      warp_out = .warpfunc_in2out(
        x = pix_grid[, 1],
        y = pix_grid[, 2],
        header_in = header_in,
        WCSref_in = WCSref_in,
        header_out = header_out,
        WCSref_out = WCSref_out,
        cores = cores
      )
    } else if (direction == 'backward') {
      warp_out = .warpfunc_out2in(
        x = pix_grid[, 1],
        y = pix_grid[, 2],
        header_in = header_in,
        WCSref_in = WCSref_in,
        header_out = header_out,
        WCSref_out = WCSref_out,
        cores = cores
      )
    }

    if(anyInfinite(warp_out)){
      warning('Infinite value in warp field, this can cause segfaults on some compilations!')
    }

    warpfield = imager::imappend(list(imager::as.cimg(matrix(
      warp_out[, 1], dim(image_out$imDat)[1], dim(image_out$imDat)[2]
    )),
    imager::as.cimg(matrix(
      warp_out[, 2], dim(image_out$imDat)[1], dim(image_out$imDat)[2]
    ))), 'c')

    rm(pix_grid)
    rm(warp_out)
  }

  image_out$imDat = imager::warp(
    im = imager::as.cimg(image_out$imDat),
    warpfield = warpfield,
    mode = switch(direction, backward = 0L, forward =
                    2L),
    interpolation = switch(
      interpolation,
      nearest = 0L,
      linear = 1L,
      cubic = 2L
    ),
    boundary_conditions = switch(
      boundary,
      dirichlet = 0L,
      neumann = 1L,
      periodic = 2L
    )
  )

  if (dofinenorm) {
    norm = matrix(1, dim(image_out$imDat)[1], dim(image_out$imDat)[2])
    norm = imager::warp(
      im = imager::as.cimg(norm),
      warpfield = warpfield,
      mode = switch(direction, backward = 0L, forward = 2L),
      interpolation = switch(
        interpolation,
        nearest = 0L,
        linear = 1L,
        cubic = 2L
      ),
      boundary_conditions = switch(
        boundary,
        dirichlet = 0L,
        neumann = 1L,
        periodic = 2L
      )
    )

    image_out$imDat = image_out$imDat / norm
    rm(norm)
  }

  if (doscale) {
    image_out$imDat = image_out$imDat * (pixscale_out / pixscale_in) ^ 2
  }

  image_out$imDat = as.matrix(image_out$imDat)

  if(dotightcrop==FALSE | keepcrop==FALSE){
    image_out = image_out[c(1L - (min_x_in - 1L), dim_out[1] - (min_x_in - 1L)),c(1L - (min_y_in - 1L), dim_out[2] - (min_y_in - 1L)), box=1] #box=1 just in case we have a single pixel left
    image_out$keyvalues$XCUTLO = 1L
    image_out$keyvalues$XCUTHI = dim_out[1]
    image_out$keyvalues$YCUTLO = 1L
    image_out$keyvalues$YCUTHI = dim_out[2]

    image_out$hdr = Rfits_keyvalues_to_hdr(image_out$keyvalues)
    image_out$header = Rfits_keyvalues_to_header(image_out$keyvalues)
    image_out$raw = Rfits_header_to_raw(Rfits_keyvalues_to_header(image_out$keyvalues))

    image_out$keynames = names(image_out$keyvalues)

    image_out$keycomments$XCUTLO = 'Low image x range'
    image_out$keycomments$XCUTHI = 'High image x range'
    image_out$keycomments$YCUTLO = 'Low image y range'
    image_out$keycomments$YCUTHI = 'High image y range'
  }else{

    if(max_x_in > dim_out[1]){
      trim_x = max_x_in - dim_out[1]
      if(dim(image_out)[1] > trim_x){
        image_out = image_out[1:(dim(image_out)[1] - trim_x), , box=1] #box=1 just in case we have a single pixel left
      }else{
        image_out = image_out[1, , box=1] #box=1 just in case we have a single pixel left
      }
      max_x_in = dim_out[1]
    }

    if(max_y_in > dim_out[2]){
      trim_y = max_y_in - dim_out[2]
      if(dim(image_out)[2] > trim_y){
        image_out = image_out[, 1:(dim(image_out)[2] - trim_y), box=1] #box=1 just in case we have a single pixel left
      }else{
        image_out = image_out[, 1, box=1] #box=1 just in case we have a single pixel left
      }
      max_y_in = dim_out[2]
    }

    if(extratight){
      final_pix = which(!is.na(image_out$imDat), arr.ind = TRUE)
      crop_x_lo = min(final_pix[,1])
      crop_x_hi = max(final_pix[,1])
      crop_y_lo = min(final_pix[,2])
      crop_y_hi = max(final_pix[,2])

      #final extra tight crop
      image_out = image_out[c(crop_x_lo, crop_x_hi), c(crop_y_lo, crop_y_hi)]

      #update where we are in the parent frame
      min_x_in = min_x_in + crop_x_lo - 1L
      max_x_in = min_x_in + dim(image_out)[1] -1L
      min_y_in = min_y_in + crop_y_lo - 1L
      max_y_in = min_y_in + dim(image_out)[2] -1L
    }

    image_out$keyvalues$XCUTLO = min_x_in
    image_out$keyvalues$XCUTHI = max_x_in
    image_out$keyvalues$YCUTLO = min_y_in
    image_out$keyvalues$YCUTHI = max_y_in

    image_out$keycomments$XCUTLO = 'Low image x range'
    image_out$keycomments$XCUTHI = 'High image x range'
    image_out$keycomments$YCUTLO = 'Low image y range'
    image_out$keycomments$YCUTHI = 'High image y range'

    image_out$keynames = names(image_out$keyvalues)

    image_out$crop = c(xlo=min_x_in, xhi=max_x_in, ylo=min_y_in, yhi=max_y_in) #we want to keep the subset location for potential later writing
  }

  if(plot){
    plot(image_out, ...)
  }

  if(warpfield_return){
    image_out$warpfield = warpfield
  }

  return(invisible(image_out))
}

propaneRebin = function(image, scale = 1,interpolation = 6){
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("The imager package is needed for this function to work. Please install it from CRAN.",
         call. = FALSE)
  }

  imdim = dim(image)
  if(scale > 1){
    scale = floor(scale)
    size_x = imdim[1]*scale - (scale - 1L)
    size_y = imdim[2]*scale - (scale - 1L)
  }else{
    scale = 1/ceiling(1/scale)
    size_x = floor(imdim[1]*scale)
    size_y = floor(imdim[2]*scale)
  }

  if(inherits(image,'Rfits_image')){
    image_resize = as.matrix(imager::resize(im=imager::as.cimg(image$imDat), size_x=size_x, size_y=size_y, interpolation_type=interpolation))
    norm = matrix(1, dim(image$imDat)[1], dim(image$imDat)[2])
    norm_resize = as.matrix(imager::resize(im=imager::as.cimg(norm), size_x=size_x, size_y=size_y, interpolation_type=interpolation))
    image_resize = (image_resize / norm_resize) / scale^2

    keyvalues_out = image$keyvalues
    keyvalues_out$NAXIS1 = dim(image_resize)[1]
    keyvalues_out$NAXIS2 = dim(image_resize)[2]
    if(scale > 1){
      keyvalues_out$CRPIX1 = (keyvalues_out$CRPIX1 - 0.5) * scale
      keyvalues_out$CRPIX2 = (keyvalues_out$CRPIX2 - 0.5) * scale
    }else{
      keyvalues_out$CRPIX1 = (keyvalues_out$CRPIX1 + 0.5) * scale
      keyvalues_out$CRPIX2 = (keyvalues_out$CRPIX2 + 0.5) * scale
    }
    keyvalues_out$CD1_1 = keyvalues_out$CD1_1 / scale
    keyvalues_out$CD1_2 = keyvalues_out$CD1_2 / scale
    keyvalues_out$CD2_1 = keyvalues_out$CD2_1 / scale
    keyvalues_out$CD2_2 = keyvalues_out$CD2_2 / scale

    image_out = list(
      imDat = image_resize,
      keyvalues = keyvalues_out,
      hdr = Rfits_keyvalues_to_hdr(keyvalues_out),
      header = Rfits_keyvalues_to_header(keyvalues_out),
      raw = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_out)),
      keynames = names(keyvalues_out),
      keycomments = as.list(rep('', length(keyvalues_out)))
    )
    names(image_out$keycomments) = image_out$keynames
    class(image_out) = c('Rfits_image', class(image_out))

    return(image_out)
  }else{
    image_resize = as.matrix(imager::resize(im=imager::as.cimg(image), size_x=size_x, size_y=size_y, interpolation_type=interpolation))
    norm = matrix(1, dim(image)[1], dim(image)[2])
    norm_resize = as.matrix(imager::resize(im=imager::as.cimg(norm), size_x=size_x, size_y=size_y, interpolation_type=interpolation))
    image_resize = (image_resize / norm_resize) / scale^2

    return(image_resize)
  }
}

propaneWarpProPane = function(propane_in, keyvalues_out=NULL, dim_out = NULL, magzero_out = NULL, ...){

  if(!is.null(magzero_out)){
    zero_point_scale = 10^(-0.4*(propane_in$image$keyvalues$MAGZERO - magzero_out))
  }else{
    magzero_out = propane_in$image$keyvalues$MAGZERO
    zero_point_scale = 1
  }

  if(!is.null(propane_in$image)){
    message('warping image')

    image_warp = propaneWarp(
      image_in = propane_in$image*zero_point_scale,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = TRUE,
      ...
    )

    image_warp$keyvalues$EXTNAME = 'image'
    image_warp$keyvalues$MAGZERO = magzero_out
  }else{
    image_warp = NULL
  }

  if(!is.null(propane_in$weight)){
    message('warping weight')

    weight_warp = propaneWarp(
      image_in = propane_in$weight,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = FALSE,
      ...
    )

    weight_warp$keyvalues$EXTNAME = 'weight'
  }else{
    weight_warp = NULL
  }

  if(!is.null(propane_in$inVar)){
    message('warping inVar')

    inVar_warp = propaneWarp(
      image_in = propane_in$inVar/(zero_point_scale^2),
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = FALSE,
      ...
    )*(Rwcs_pixscale(propane_in$inVar$keyvalues)^4 / Rwcs_pixscale(keyvalues_out)^4)

    inVar_warp$keyvalues$EXTNAME = 'inVar'
    inVar_warp$keyvalues$MAGZERO = magzero_out
  }else{
    inVar_warp = NULL
  }

  if(!is.null(propane_in$exp)){
    message('warping exp')

    exp_warp = propaneWarp(
      image_in = propane_in$exp,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = FALSE,
      ...
    )

    exp_warp$keyvalues$EXTNAME = 'exp'
  }else{
    exp_warp = NULL
  }

  if(!is.null(propane_in$cold)){
    message('warping cold')

    cold_warp = propaneWarp(
      image_in = propane_in$cold*zero_point_scale,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = TRUE,
      ...
    )

    cold_warp$keyvalues$EXTNAME = 'cold'
    cold_warp$keyvalues$MAGZERO = magzero_out
  }else{
    cold_warp = NULL
  }

  if(!is.null(propane_in$hot)){
    message('warping hot')

    hot_warp = propaneWarp(
      image_in = propane_in$hot*zero_point_scale,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = TRUE,
      ...
    )

    hot_warp$keyvalues$EXTNAME = 'hot'
    hot_warp$keyvalues$MAGZERO = magzero_out
  }else{
    hot_warp = NULL
  }

  if(!is.null(propane_in$clip)){
    message('warping clip')

    clip_warp = propaneWarp(
      image_in = propane_in$clip,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = FALSE,
      ...
    )

    clip_warp$keyvalues$EXTNAME = 'clip'
  }else{
    clip_warp = NULL
  }

  output = list(
    image = image_warp,
    weight = weight_warp,
    inVar = inVar_warp,
    exp = exp_warp,
    cold = cold_warp,
    hot = hot_warp,
    clip=clip_warp
  )

  class(output) = "ProPane"
  return(invisible(output))
}
