propaneBadPix = function(image, mask=NULL, inVar=NULL, smooth=1, sigma=10, pixcut=1, cold=FALSE, hot=TRUE,
                         dilate=FALSE, size=3, return='image', patch=FALSE, allow_write=FALSE, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  assert(checkClass(image,'Rfits_image'), checkClass(image,'Rfits_pointer'), checkClass(image,'matrix'))
  assertMatrix(mask, null.ok=TRUE)
  assertMatrix(inVar, null.ok=TRUE)
  assertNumeric(smooth, len=1)
  assertNumeric(sigma, max.len=2)
  assertIntegerish(pixcut, max.len=2)
  assertFlag(cold)
  assertFlag(hot)
  assertLogical(dilate, max.len=2)
  assertIntegerish(size, max.len=2)
  assertChoice(return, c('image', 'mask', 'loc'))
  assertFlag(patch)
  assertFlag(allow_write)
  assertFlag(plot)

  if(length(sigma) == 1){
    sigma = rep(sigma, 2)
  }

  if(length(pixcut) == 1){
    pixcut = rep(pixcut, 2)
  }

  if(length(dilate) == 1){
    dilate = rep(dilate, 2)
  }

  if(length(size) == 1){
    size = rep(size, 2)
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

  image_blur = as.matrix(imager::isoblur(imager::as.cimg(image_data),smooth,na.rm=TRUE))
  image_diff = image_data - image_blur

  if(is.null(inVar)){
    quancuts = quantile(image_diff, c(0.1586553, 0.5), na.rm=TRUE)
  }else{
    quancuts = c(0,1)
    image_diff = image_diff*sqrt(inVar)
  }

  image_diff[!is.finite(image_diff)] = NA

  if(cold){
    thresh_cold = quancuts[2] - (quancuts[2] - quancuts[1])*sigma[1]
  }else{
    thresh_cold = -Inf
  }

  if(hot){
    thresh_hot = quancuts[2] + (quancuts[2] - quancuts[1])*sigma[2]
  }else{
    thresh_hot = Inf
  }

  if(!is.null(mask)){
    image_diff[mask_in_sel] = 0L
    image_data[mask_in_sel] = image_orig[mask_in_sel]
  }

  mask_rough_cold = (image_diff < thresh_cold) & (image_data <= 0)
  mask_label_cold = as.matrix(imager::label(imager::as.cimg(mask_rough_cold)))
  sel_cold = which(tabulate(mask_label_cold) <= pixcut[1])
  mask_rough_cold[!mask_label_cold %in% sel_cold] = 0L

  mask_rough_hot = (image_diff > thresh_hot) & (image_data >= 0)
  mask_label_hot = as.matrix(imager::label(imager::as.cimg(mask_rough_hot)))
  sel_hot = which(tabulate(mask_label_hot) <= pixcut[2])
  mask_rough_hot[!mask_label_hot %in% sel_hot] = 0L

  if(dilate[1]){
    if(size[1] %% 2 == 0){
      size[1] = size[1] + 1
    }
    kernel = imager::px.circle(size[1]/2, size[1], size[1])
    mask_rough_cold = as.matrix(imager::dilate(imager::as.cimg(mask_rough_cold), mask=kernel))
  }

  if(dilate[2]){
    if(size[2] %% 2 == 0){
      size[2] = size[2] + 1
    }
    kernel = imager::px.circle(size[2]/2, size[2], size[2])
    mask_rough_hot = as.matrix(imager::dilate(imager::as.cimg(mask_rough_hot), mask=kernel))
  }

  mask_loc = which(mask_rough_cold + mask_rough_hot > 0L, arr.ind=TRUE)

  if(return == 'image'){
    image_data[mask_loc] = NA

    if(patch){
      image_blur = as.matrix(imager::isoblur(imager::as.cimg(image_data),smooth,na.rm=TRUE))
      image_data[mask_loc] = image_blur[mask_loc]
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

propanePatchPix = function(image, mask=NULL, smooth=1, dilate=FALSE, size=3, allow_write=FALSE,
                           plot=FALSE, ...){
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

  if(dilate){
    if(size %% 2 == 0){
      size = size + 1
    }
    kernel = imager::px.circle(size/2,size,size)
    mask_rough = as.matrix(imager::dilate(imager::as.cimg(is.na(image_data)), mask=kernel))
    image_data[mask_rough > 0L] = NA
  }

  image_blur = as.matrix(imager::isoblur(imager::as.cimg(image_data),smooth,na.rm=TRUE))

  if(!is.null(mask)){
    image_data[mask_in_sel] = image_orig[mask_in_sel]
    sel_patch = which(is.na(image_data) & !mask_in_sel, arr.ind=TRUE)
  }else{
    sel_patch = which(is.na(image_data), arr.ind=TRUE)
  }

  image_data[sel_patch] = image_blur[sel_patch]

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
