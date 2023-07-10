propaneMaskPix = function(image, smooth=1, sigma=10, pixcut=1, cold=FALSE, hot=TRUE, return='image',
                         patch=FALSE, allow_write=FALSE, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  if(inherits(image, 'Rfits_image')){
    image_data = image$imDat
  }else if(inherits(image, 'Rfits_pointer')){
    if(allow_write){
      image_data = image[,,header=FALSE]
    }else{
      image = image[,]
      image_data = image$imDat
    }
  }else{
    image_data = image
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
    }else{
      image = image_data
    }

    if(plot){
      plot(image, ...)
    }
    return(invisible(image))
  }else if(return == 'mask'){
    mask = matrix(0L, dim(image)[1], dim(image)[2])
    mask[mask_loc] = 1L
    if(plot){
      magimage(mask, ...)
    }
    return(invisible(mask))
  }else if(return == 'loc'){
    return(invisible(mask_loc))
  }
}
