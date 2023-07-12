propaneBadPix = function(image, mask=NULL, smooth=1, sigma=10, pixcut=1, cold=FALSE, hot=TRUE, return='image',
                         patch=FALSE, allow_write=FALSE, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  assert(checkClass(image,'Rfits_image'), checkClass(image,'Rfits_pointer'), checkClass(image,'matrix'))
  assertNumeric(smooth, len=1)
  assertNumeric(sigma, len=1)
  assertIntegerish(pixcut, len=1)
  assertFlag(cold)
  assertFlag(hot)
  assertChoice(return, c('image', 'mask', 'loc'))
  assertFlag(patch)
  assertFlag(allow_write)
  assertFlag(plot)

  if(inherits(image, 'Rfits_image')){
    image_data = image$imDat
  }else if(inherits(image, 'Rfits_pointer')){
    if(allow_write){
      image_data = image[,,header=FALSE]
    }else{
      image = image[,]
      image_data = image$imDat
    }
  }else if(is.matrix(image)){
    image_data = image
  }else{
    stop('image input must be Rfits_image, Rfits_pointer of matrix!')
  }

  if(!is.null(mask)){
    mask_in_sel = which(mask > 0L)
    image_orig = image_data
    image_data[mask_in_sel] = NA
  }

  blur = as.matrix(imager::isoblur(imager::as.cimg(image_data),smooth,na.rm=TRUE))
  image_diff = image_data - blur

  quancuts = quantile(image_diff, c(0.1586553, 0.5), na.rm=TRUE)

  if(cold){
    thresh_cold = quancuts[2] - (quancuts[2] - quancuts[1])*sigma
  }else{
    thresh_cold = -Inf
  }

  if(hot){
    thresh_hot = quancuts[2] + (quancuts[2] - quancuts[1])*sigma
  }else{
    thresh_hot = Inf
  }

  if(!is.null(mask)){
    image_diff[mask_in_sel] = 0L
    image_data[mask_in_sel] = image_orig[mask_in_sel]
  }

  mask_rough = (image_diff < thresh_cold) | (image_diff > thresh_hot)
  mask_label = as.matrix(imager::label(imager::as.cimg(mask_rough)))
  sel = which(tabulate(mask_label) <= pixcut)
  mask_loc = which(matrix(mask_label %in% sel, dim(image)[1], dim(image)[2]), arr.ind=TRUE)

  if(return == 'image'){
    image_data[mask_loc] = NA

    if(patch){
      blur = as.matrix(imager::isoblur(imager::as.cimg(image_data),smooth,na.rm=TRUE))
      image_data[mask_loc] = blur[mask_loc]
    }

    if(inherits(image, 'Rfits_image')){
      image$imDat = image_data
    }else if(inherits(image, 'Rfits_pointer')){
      image[mask_loc, allow_write=allow_write] = image_data[mask_loc]
    }else if(is.matrix(image)){
      image = image_data
    }

    if(plot){
      if(is.matrix(image)){
        magimage(image, ...)
      }else{
        plot(image, ...)
      }
    }

    return(invisible(image))
  }else if(return == 'mask'){
    mask_out = matrix(0L, dim(image)[1], dim(image)[2])
    mask_out[mask_loc] = 1L
    if(plot){
      magimage(mask_out, ...)
    }
    return(invisible(mask_out))
  }else if(return == 'loc'){
    return(invisible(mask_loc))
  }
}

propanePatchPix = function(image, mask=NULL, smooth=1, allow_write=FALSE, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  assert(checkClass(image,'Rfits_image'), checkClass(image,'Rfits_pointer'), checkClass(image,'matrix'))

  if(inherits(image, 'Rfits_image')){
    image_data = image$imDat
  }else if(inherits(image, 'Rfits_pointer')){
    if(allow_write){
      image_data = image[,,header=FALSE]
    }else{
      image = image[,]
      image_data = image$imDat
    }
  }else if(is.matrix(image)){
    image_data = image
  }else{
    stop('image input must be Rfits_image, Rfits_pointer of matrix!')
  }

  if(!is.null(mask)){
    mask_in_sel = mask > 0L
    image_orig = image_data
    image_data[mask_in_sel] = NA
  }

  blur = as.matrix(imager::isoblur(imager::as.cimg(image_data),smooth,na.rm=TRUE))

  if(!is.null(mask)){
    image_data[mask_in_sel] = image_orig[mask_in_sel]
  }

  sel_patch = which(is.na(image_data) & !mask_in_sel, arr.ind=TRUE)
  image_data[sel_patch] = blur[sel_patch]

  if(inherits(image, 'Rfits_image')){
    image$imDat = image_data
  }else if(inherits(image, 'Rfits_pointer')){
    image[sel_patch, allow_write=allow_write] = image_data[sel_patch]
  }else if(is.matrix(image)){
    image = image_data
  }

  if(plot){
    if(is.matrix(image)){
      magimage(image, ...)
    }else{
      plot(image, ...)
    }
  }
}
