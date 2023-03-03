.get_pix = function(image){
  if(is.matrix(image)){
    return(image)
  }else if(inherits(image,'Rfits_image')){
    return(image$imDat)
  }else if(inherits(image,'Rfits_pointer')){
    return(image[,]$imDat)
  }else{
    stop('Must be of class matrix, Rfits_image, Rfits_pointer')
  }
}

propanePatch = function(image_inVar, image_med, diff_type='scale', threshold=5,
                        scale_type='quan', scale_val=0.3, hot=TRUE, cold=TRUE){
  if(!is.null(image_inVar$keyvalues)){
    keyvalues = image_inVar$keyvalues
  }else{
    keyvalues = NULL
  }
  image_inVar = .get_pix(image_inVar)
  image_med = .get_pix(image_med)
  if(hot & cold){
    im_diff = abs(image_inVar - image_med)
  }else if(hot){
    im_diff = image_inVar - image_med
  }else if (cold){
    im_diff = image_med - image_inVar
  }else{
    stop('Doing nothing! Need to clip one or both of hot/cold')
  }

  if(diff_type=='scale'){
    if(scale_type == 'quan'){
      im_diff = im_diff/quantile(abs(image_med), probs=scale_val, na.rm=TRUE)
    }else if(scale_type == 'median' | scale_type == 'med'){
      im_diff = im_diff/median(abs(image_med), na.rm=TRUE)
    }else if(scale_type == 'num' | scale_type == 'val'){
      im_diff = im_diff/scale_val
    }
  }else if(diff_type=='rel'){
    if(scale_type == 'quan'){
      ref_val = quantile(abs(image_med), probs=scale_val, na.rm=TRUE)
    }else if(scale_type == 'median' | scale_type == 'med'){
      ref_val = median(abs(image_med), na.rm=TRUE)
    }else if(scale_type == 'num' | scale_type == 'val'){
      ref_val = scale_val
    }
    denom = abs(image_med)
    denom[denom < ref_val] = ref_val
    im_diff = im_diff / denom
  }else if(diff_type=='diff'){
    #Doing nothing
  }else{
    stop('Doing nothing! type must be one of diff, scale, rel')
  }

  sel_threshold = (im_diff > threshold)
  sel_NA = is.na(image_inVar)
  sel_both = which(sel_threshold | sel_NA)
  output = image_inVar
  output[sel_both] = image_med[sel_both]

  if(!is.null(keyvalues)){
    output = Rfits_create_image(output, keyvalues=keyvalues)
  }

  message('Patching ',
          round(100*length(which(sel_threshold))/prod(dim(image_inVar)), digits=2),'% of threshold pixels and ',
          round(100*length(which(sel_NA))/prod(dim(image_inVar)), digits=2),'% of NA pixels'
          )

  patch_mat = matrix(0L, dim(image_inVar)[1], dim(image_inVar)[2])
  patch_mat[which(sel_NA)] = 1L
  patch_mat[which(sel_threshold)] = 2L
  if(!is.null(keyvalues)){
    keyvalues$EXTNAME = 'patch'
    patch_mat = Rfits_create_image(patch_mat, keyvalues=keyvalues)
  }
  return(list(image=output, patch=patch_mat))
}
