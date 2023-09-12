propaneTweak = function(image_ref, image_pre_fix, delta_max=c(3,0), Nmeta=3, quan_cut=c(0.98,0.9999),
                        stretch='lin', WCS_match=TRUE, cores=1, shift_int=TRUE, algotype='optim',
                        Niter=1e4, return_image=TRUE, direction='backward', final_centre=FALSE,
                        cutcheck=FALSE, quick=FALSE, verbose=TRUE){

  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
  }

  timestart = proc.time()[3]

  algotype = tolower(algotype)

  if(shift_int){
    WCS_match = TRUE
  }

  if(WCS_match){
    if(all(dim(image_ref) != dim(image_pre_fix))){
      stop('If WCS_match or shift_int is TRUE then dimensions of image_ref and image_pre_fix must match!')
    }
  }

  if(WCS_match==FALSE){
    if(inherits(image_ref, 'Rfits_image') == FALSE){
      stop('image_ref must be class Rfits_image if WCS_match=FALSE')
    }
    if(inherits(image_pre_fix, 'Rfits_image') == FALSE){
      stop('image_pre_fix must be class Rfits_image if WCS_match=FALSE')
    }
  }

  if(inherits(image_ref, 'Rfits_image')){
    if(WCS_match){
      image_ref = image_ref$imDat
    }
  }

  image_pre_fix_orig = image_pre_fix

  if(inherits(image_pre_fix, 'Rfits_image')){
    if(WCS_match){
      image_pre_fix = image_pre_fix$imDat
    }
  }

  #dim_orig = dim(image_ref)

  if(verbose){
    message("Trimming comparison region and auto scaling")
  }

  if(WCS_match){
    trim_lim = which((!is.na(image_ref)) & (!is.na(image_pre_fix)), arr.ind = TRUE)
    x_lo = min(trim_lim[,1])
    x_hi = max(trim_lim[,1])
    y_lo = min(trim_lim[,2])
    y_hi = max(trim_lim[,2])

    image_ref = image_ref[x_lo:x_hi, y_lo:y_hi]
    image_pre_fix = image_pre_fix[x_lo:x_hi, y_lo:y_hi]
  }else{
    #this might look backwards, but we don't want to warp image_pre_fix here, just maximally crop it.
    image_overlap = propaneCropOverlap(image_ref=image_pre_fix, image_warp=image_ref)
    image_pre_fix = image_overlap$image_ref
    image_ref_warp = image_overlap$image_warp
  }

  if(inherits(image_pre_fix, 'Rfits_image')){
    if(WCS_match){
      im_med = median(image_pre_fix$imDat, na.rm=TRUE)
      im_quan_lo = quantile(image_pre_fix$imDat, quan_cut[1], na.rm=TRUE)
    }else{
      sel = which(!is.na(image_ref_warp$imDat))
      im_med = median(image_pre_fix$imDat[sel], na.rm=TRUE)
      im_quan_lo = quantile(image_pre_fix$imDat[sel], quan_cut[1], na.rm=TRUE)
    }
    if(length(quan_cut)==2){
      image_pre_fix$imDat[image_pre_fix$imDat > quantile(image_pre_fix$imDat, quan_cut[2], na.rm=TRUE)] = 0
    }
    image_pre_fix$imDat = image_pre_fix$imDat - im_med
    image_pre_fix$imDat = image_pre_fix$imDat / im_quan_lo
    image_pre_fix$imDat[image_pre_fix$imDat < 1] = 0
    image_pre_fix$imDat = .stretch(image_pre_fix$imDat, stretch=stretch)
    if(cutcheck){
      plot(image_pre_fix)
      legend('topleft', 'image_pre_fix')
    }
  }else{
    if(length(quan_cut)==2){
      image_pre_fix[image_pre_fix > quantile(image_pre_fix, quan_cut[2], na.rm=TRUE)] = 0
    }
    image_pre_fix = image_pre_fix - median(image_pre_fix, na.rm=TRUE)
    image_pre_fix = image_pre_fix / quantile(image_pre_fix, quan_cut[1], na.rm=TRUE)
    image_pre_fix[image_pre_fix < 1] = 0
    image_pre_fix = .stretch(image_pre_fix, stretch=stretch)
    if(cutcheck){
      magimage(image_pre_fix)
      legend('topleft', 'image_pre_fix')
    }
  }

  if(inherits(image_ref, 'Rfits_image')){
    if(length(quan_cut)==2){
      image_ref$imDat[image_ref$imDat > quantile(image_ref$imDat, quan_cut[2], na.rm=TRUE)] = NA
    }
    image_ref$imDat = image_ref$imDat - median(image_ref$imDat, na.rm=TRUE)
    image_ref$imDat = image_ref$imDat / quantile(image_ref$imDat, quan_cut[1], na.rm=TRUE)
    image_ref$imDat[image_ref$imDat < 1] = NA
    image_ref$imDat = .stretch(image_ref$imDat, stretch=stretch)
    if(cutcheck){
      plot(image_ref)
      legend('topleft', 'image_ref')
      return(NULL)
    }
  }else{
    if(length(quan_cut)==2){
      image_ref[image_ref > quantile(image_ref, quan_cut[2], na.rm=TRUE)] = NA
    }
    image_ref = image_ref - median(image_ref, na.rm=TRUE)
    image_ref = image_ref / quantile(image_ref, quan_cut[1], na.rm=TRUE)
    image_ref[image_ref < 1] = NA
    image_ref = .stretch(image_ref, stretch=stretch)
    if(cutcheck){
      magimage(image_ref)
      legend('topleft', 'image_ref')
      return(NULL)
    }
  }


  #scale = 1 # not really used anymore
  i = NULL #to avoid warnings

  #pix_cost_use = which(!is.na(image_ref) & !is.na(image_pre_fix))

  if(verbose){
    message('Optimising tweak')
  }

  if(shift_int){

    registerDoParallel(cores=cores)

    Ndelta = length(-delta_max[1]:delta_max[1])

    for(Nshift in 1:Nmeta){

      if(verbose){
        message('  Meta shift ', Nshift, ' of ', Nmeta)
      }

      if(Nshift == 1){
        current_par = c(0L, 0L)
      }

      grid_search = expand.grid(-delta_max[1]:delta_max[1] + current_par[1], -delta_max[1]:delta_max[1] + current_par[2])

      cost_mat = foreach(i = 1:dim(grid_search)[1], .combine='c')%dopar%{
        cost = .mat_diff_sum(image_ref, image_pre_fix, 1, grid_search[i,1], grid_search[i,2]) #1 is the unused scale param
        if(verbose){
          message(grid_search[i,1],' ', grid_search[i,2],' ',cost)
        }
        return(cost)
      }

      current_par = as.integer(grid_search[which.min(cost_mat),])

      if(verbose){
        message('    Current best: ', current_par[1], ' ', current_par[2],': ',round(min(cost_mat)))
      }

      at_lim = current_par[1] == min(grid_search[,1]) | current_par[1] == max(grid_search[,1]) | current_par[2] == min(grid_search[,2]) | current_par[2] == max(grid_search[,2])
      if(at_lim == FALSE){
        if(final_centre){
          if(verbose){
            message('  Final centre shift')
          }

          grid_search = expand.grid(-delta_max[1]:delta_max[1] + current_par[1], -delta_max[1]:delta_max[1] + current_par[2])

          cost_mat = foreach(i = 1:dim(grid_search)[1], .combine='c')%dopar%{
            cost = .mat_diff_sum(image_ref, image_pre_fix, 1, grid_search[i,1], grid_search[i,2]) #1 is the unused scale param
            if(verbose){
              message(grid_search[i,1],' ', grid_search[i,2],' ',cost)
            }
            return(cost)
          }

          current_par = as.integer(grid_search[which.min(cost_mat),])

          break
        }else{
          break
        }
      }else{
        if(verbose){
          message('    Limit hit!')
        }
      }
    }

    optim_out = list()
    optim_out$par = current_par
    optim_out$value = min(cost_mat)
    optim_out$counts = Ndelta^2

    cos_weight = exp(-(cost_mat - min(cost_mat))/2)
    par_weight_x = sum(grid_search[,1] * cos_weight) / sum(cos_weight)
    par_weight_y = sum(grid_search[,2] * cos_weight) / sum(cos_weight)
    optim_out$par_weight = c(par_weight_x, par_weight_y)

  }else{

    if(delta_max[2] == 0){
      par = c(0,0)
      lower = c(-delta_max[1], -delta_max[1])
      upper = c(delta_max[1], delta_max[1])
    }else{
      par = c(0,0,0)
      lower = c(-delta_max[1], -delta_max[1], -delta_max[2])
      upper = c(delta_max[1], delta_max[1], delta_max[2])
    }

    if(algotype == 'cma'){
      if(!requireNamespace("cmaes", quietly = TRUE)){
        stop('The cmaes package is needed for matching to work. Please install from CRAN.', call. = FALSE)
      }
    }

    if(quick==FALSE){
      if(algotype == 'optim'){
        optim_out = optim(par = par,
                          fn = .cost_fn_image,
                          method = "L-BFGS-B",
                          lower = lower,
                          upper = upper,
                          image_ref = image_ref,
                          image_pre_fix = image_pre_fix,
                          direction = direction,
                          pix_cost_use = NULL,
                          shift_int = FALSE,
                          WCS_match = WCS_match,
                          cores = cores
        )
      }else if(algotype == 'cma'){
        optim_out = cmaes::cma_es(par = par,
                          fn = .cost_fn_image,
                          method = "L-BFGS-B",
                          lower = lower,
                          upper = upper,
                          image_ref = image_ref,
                          image_pre_fix = image_pre_fix,
                          direction = direction,
                          pix_cost_use = NULL,
                          shift_int = FALSE,
                          WCS_match = WCS_match,
                          control = list(maxit = Niter),
                          cores = cores
        )
      }
    }else{
      if(inherits(image_pre_fix, 'Rfits_image')){
        image_pre_fix_xysub = which(image_pre_fix$imDat > 0, arr.ind = TRUE)
        image_pre_fix_xysub = cbind(image_pre_fix_xysub, image_pre_fix$imDat[image_pre_fix_xysub])
      }else{
        image_pre_fix_xysub = which(image_pre_fix > 0, arr.ind = TRUE)
        image_pre_fix_xysub = cbind(image_pre_fix_xysub, image_pre_fix[image_pre_fix_xysub])
      }

      #image_pre_fix_xysub[,1:2] = image_pre_fix_xysub[,1:2] - 0.5 #to get into R units

      if(algotype == 'optim'){
        optim_out = optim(par = par,
                          fn = .cost_fn_image_approx,
                          method = "L-BFGS-B",
                          lower = lower,
                          upper = upper,
                          image_ref = image_ref,
                          image_pre_fix_xysub = image_pre_fix_xysub,
                          keyvalues_pre_fix = image_pre_fix$keyvalues,
                          WCS_match = WCS_match,
                          xcen_rot = mean(image_pre_fix_xysub[,1], na.rm=TRUE),
                          ycen_rot = mean(image_pre_fix_xysub[,2], na.rm=TRUE),
                          cores = cores
        )
      }else if(algotype == 'cma'){
        optim_out = cmaes::cma_es(
          par = par,
          fn = .cost_fn_image_approx,
          lower = lower,
          upper = upper,
          image_ref = image_ref,
          image_pre_fix_xysub = image_pre_fix_xysub,
          keyvalues_pre_fix = image_pre_fix$keyvalues,
          WCS_match = WCS_match,
          xcen_rot = mean(image_pre_fix_xysub[,1], na.rm=TRUE),
          ycen_rot = mean(image_pre_fix_xysub[,2], na.rm=TRUE),
          control = list(maxit = Niter),
          cores = cores
        )
      }else{
        stop('')
      }
    }
  }

  if(return_image){
    if(verbose){
      message('Creating output image_post_fix')
    }

    if(inherits(image_pre_fix_orig, 'Rfits_image')){
      if(WCS_match){
        image_post_fix = image_pre_fix_orig

        if(optim_out$par[1] != 0 | optim_out$par[2] != 0){
          image_post_fix$imDat = .cost_fn_image(par = optim_out$par,
                                    image_ref = image_ref,
                                    image_pre_fix = image_pre_fix_orig$imDat,
                                    direction = direction,
                                    return = 'image_post_fix',
                                    shift_int = FALSE,
                                    WCS_match = TRUE,
                                    cores = cores)
        }
      }else{
        if(length(optim_out$par) == 2){
          image_post_fix = propaneWCSmod(image_pre_fix_orig, delta_x=optim_out$par[1], delta_y=optim_out$par[2], recen=FALSE)
        }else{
          image_post_fix = propaneWCSmod(image_pre_fix_orig, delta_x=optim_out$par[1], delta_y=optim_out$par[2], delta_rot=optim_out$par[3], recen=FALSE)
        }
      }
    }else{
      if(optim_out$par[1] != 0 | optim_out$par[2] != 0){
        image_post_fix = .cost_fn_image(par = optim_out$par,
                                       image_ref = image_ref,
                                       image_pre_fix = image_pre_fix_orig,
                                       direction = direction,
                                       return = 'image_post_fix',
                                       shift_int = FALSE,
                                       WCS_match = TRUE,
                                       cores = cores)
      }else{
        image_post_fix = image_pre_fix_orig
      }

      # image_post_fix = matrix(NA, dim_orig[1], dim_orig[2])
      # image_post_fix[x_lo:x_hi, y_lo:y_hi] = image_post_fix_temp
    }
  }else{
    image_post_fix = NULL
  }

  if(shift_int){
    return(invisible(list(optim_out = optim_out, image_post_fix = image_post_fix, time = proc.time()[3] - timestart, cost_mat=list(x= -delta_max[1]:delta_max[1] + current_par[1], y= -delta_max[1]:delta_max[1] + current_par[2], z=matrix(cost_mat,Ndelta,Ndelta)), at_lim=at_lim)))
  }else{
    return(invisible(list(optim_out = optim_out, image_post_fix = image_post_fix, time = proc.time()[3] - timestart)))
  }
}

propaneTweakImage = propaneTweak

propaneTweakCat = function(cat_ref, cat_pre_fix, delta_max=c(100,0), mode='pix',
                           keyvalues_pre_fix=NULL){

  if(mode == 'coord'){
    if(is.null(keyvalues_pre_fix)){
      stop('keyvalues_pre_fix must be provided when input mode is coord!')
    }
    cat_ref = Rwcs_s2p(cat_ref, keyvalues=keyvalues_pre_fix)
    cat_pre_fix = Rwcs_s2p(cat_ref, cat_pre_fix=keyvalues_pre_fix)
  }

  if(!requireNamespace("RANN", quietly = TRUE)){
    stop('The RANN package is needed for matching to work. Please install from CRAN.', call. = FALSE)
  }

  if(delta_max[2] == 0){
    par = c(0,0)
    lower = c(-delta_max[1], -delta_max[1])
    upper = c(delta_max[1], delta_max[1])
    xcen_rot = 0
    ycen_rot = 0
  }else{
    xcen_rot = mean(cat_ref[,1], na.rm=TRUE)
    ycen_rot = mean(cat_ref[,2], na.rm=TRUE)

    par = c(0,0,0)
    lower = c(-delta_max[1], -delta_max[1], -delta_max[2])
    upper = c(delta_max[1], delta_max[1], delta_max[2])
  }

  optim_out = optim(par = par,
                    fn = .cost_fn_cat,
                    method = "L-BFGS-B",
                    lower = lower,
                    upper = upper,
                    cat_ref = cat_ref,
                    cat_pre_fix = cat_pre_fix,
                    rad_max = delta_max[1],
                    xcen_rot = xcen_rot,
                    ycen_rot = ycen_rot
  )

  return(optim_out)
}

# propaneNifty = function(image_ref, image_pre_fix, scope='rigid', box=1900, grid=box, buffer=30, cores=4, ...){
#   if(length(box) == 1){
#     box = rep(box, 2)
#     if(missing(grid)){grid = box}
#   }
#
#   if(length(grid) == 1){
#     grid = rep(grid ,2)
#   }
#
#   if(box[1] > dim(image_ref)[1]){box[1] = dim(image_ref)[1]}
#   if(box[2] > dim(image_ref)[2]){box[2] = dim(image_ref)[2]}
#   if(grid[1] > dim(image_ref)[1]){grid[1] = dim(image_ref)[1]}
#   if(grid[2] > dim(image_ref)[2]){grid[2] = dim(image_ref)[2]}
#
#   xdim_ref = dim(image_ref)[1]
#   ydim_ref = dim(image_ref)[2]
#   xdim_pre_fix = dim(image_pre_fix)[1]
#   ydim_pre_fix = dim(image_pre_fix)[2]
#
#   tweak_grid_ref = expand.grid(seq(1L,xdim_ref,by=grid[1]), seq(1L,ydim_ref,by=grid[2]))
#   tweak_grid_ref = cbind(tweak_grid_ref[,1], tweak_grid_ref[,1] + box[1] - 1L, tweak_grid_ref[,2], tweak_grid_ref[,2] + box[2] - 1L)
#   tweak_grid_ref[tweak_grid_ref[,2] > xdim_ref, 2] = xdim_ref
#   tweak_grid_ref[tweak_grid_ref[,4] > ydim_ref, 4] = ydim_ref
#
#   tweak_grid_pre_fix = tweak_grid_ref
#   tweak_grid_pre_fix[,1] = tweak_grid_pre_fix[,1] - buffer
#   tweak_grid_pre_fix[,2] = tweak_grid_pre_fix[,2] + buffer
#   tweak_grid_pre_fix[,3] = tweak_grid_pre_fix[,3] - buffer
#   tweak_grid_pre_fix[,4] = tweak_grid_pre_fix[,4] + buffer
#   tweak_grid_pre_fix[tweak_grid_pre_fix[,1] < 1L, 1] = 1L
#   tweak_grid_pre_fix[tweak_grid_pre_fix[,2] > ydim_pre_fix, 2] = xdim_pre_fix
#   tweak_grid_pre_fix[tweak_grid_pre_fix[,3] < 1L, 3] = 1L
#   tweak_grid_pre_fix[tweak_grid_pre_fix[,4] > ydim_pre_fix, 4] = ydim_pre_fix
#
#   registerDoParallel(cores=cores)
#
#   tweak_out = foreach(i = 1:dim(tweak_grid_ref)[1])%dopar%{
#     image_pre_fix_cut = image_pre_fix[tweak_grid_pre_fix[i,1]:tweak_grid_pre_fix[i,2], tweak_grid_pre_fix[i,3]:tweak_grid_pre_fix[i,4]]
#     image_ref_cut = image_ref[tweak_grid_ref[i,1]:tweak_grid_ref[i,2], tweak_grid_ref[i,3]:tweak_grid_ref[i,4]]
#     sourceMask = !is.na(image_pre_fix_cut)
#     targetMask = !is.na(image_ref_cut)
#     local_tweak = try(niftyreg(
#           source = image_pre_fix_cut,
#           target = image_ref_cut,
#           scope = scope,
#           sourceMask = sourceMask,
#           targetMask = targetMask,
#           ...
#         ))
#     if(inherits(local_tweak, 'try-error')){
#       return(NULL)
#     }else{
#       #local_tweak$image_post_fix_cut = local_tweak$image
#       #local_tweak$image_pre_fix_cut = image_pre_fix_cut
#       #local_tweak$image_ref_cut = image_ref_cut
#       return(local_tweak)
#     }
#   }
#
#   image_post_fix = matrix(NA, xdim_ref, ydim_ref)
#
#   for(i in 1:dim(tweak_grid_ref)[1]){
#     if(!is.null(tweak_out[[i]])){
#       image_post_fix[tweak_grid_ref[i,1]:tweak_grid_ref[i,2], tweak_grid_ref[i,3]:tweak_grid_ref[i,4]] = tweak_out[[i]]$image
#     }
#   }
#
#   return(list(image_post_fix=image_post_fix, tweak_grid_ref=tweak_grid_ref, tweak_grid_pre_fix=tweak_grid_pre_fix, tweak_out=tweak_out))
# }

propaneTran = function(image, delta_x = 0, delta_y = 0, delta_rot = 0, xcen_rot = dim(image)[1]/2 + 0.5,
                       ycen_rot = dim(image)[1]/2 + 0.5, direction='backward', padNA=TRUE, shift_int=TRUE){
  #delta refers to the direction we shift the image, not the view point.
  #postive delta_x moves image on our viewer to the right
  #positive delta_y moves image on our viewer up
  #positive delta_rot rotates image clockwise on our viewer

  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
  }

  if(delta_rot != 0){
    shift_int = FALSE
  }

  if(shift_int){
    if(padNA){
      min_val = min(image, na.rm=TRUE)
      image = image - min_val + 1
      image = as.matrix(imager::imshift(imager::as.cimg(image), delta_x=round(delta_x), delta_y=round(delta_y)))
      image[image == 0] = NA
      image = image + min_val - 1
      return(image)
    }else{
      return(as.matrix(imager::imshift(imager::as.cimg(image), delta_x=round(delta_x), delta_y=round(delta_y))))
    }
  }else{
    if(direction == 'forward'){
      formals(.map.tran)$delta_x = delta_x
      formals(.map.tran)$delta_y = delta_y
      formals(.map.tran)$delta_rot = delta_rot
      formals(.map.tran)$xcen_rot = xcen_rot
      formals(.map.tran)$ycen_rot = ycen_rot
    }
    if(direction == 'backward'){
      formals(.map.tran)$delta_x = -delta_x
      formals(.map.tran)$delta_y = -delta_y
      formals(.map.tran)$delta_rot = -delta_rot
      formals(.map.tran)$xcen_rot = xcen_rot
      formals(.map.tran)$ycen_rot = ycen_rot
    }
    local_fun = function(x, y){
      .map.tran(x = x, y = y)
    }
    norm_mat = matrix(1, dim(image)[1], dim(image)[2])
    norm_mat = as.matrix(imager::imwarp(imager::as.cimg(norm_mat), map=local_fun, direction=direction, coordinates='absolute'))
    if(padNA){
      min_val = min(image, na.rm=TRUE)
      image = image - min_val + 1
      image = as.matrix(imager::imwarp(imager::as.cimg(image), map=local_fun, direction=direction, coordinates='absolute'))
      image[image == 0] = NA
      image = image/norm_mat
      image = image + min_val - 1
      return(image)
    }else{
      return(as.matrix(imager::imwarp(imager::as.cimg(image), map=local_fun, direction=direction, coordinates='absolute'))/norm_mat)
    }
  }
}

propaneWCSmod = function(input, delta_x = 0, delta_y = 0, delta_rot = 0, recen = FALSE){

  if(inherits(input, 'Rfits_keylist')){
    keyvalues = input
    keylist = TRUE
  }else if(inherits(input, c('Rfits_image', 'Rfits_pointer', 'Rfits_header'))){
    keyvalues = input$keyvalues
    keylist = FALSE
  }else{
    stop('input must be one of Rfits_keylist, Rfits_image, Rfits_pointer, Rfits_header')
  }

  if(delta_x != 0){
    keyvalues$CRPIX1 = keyvalues$CRPIX1 - delta_x
  }

  if(delta_y != 0){
    keyvalues$CRPIX2 = keyvalues$CRPIX2 - delta_y
  }

  # I think the above is the more accurate way to do it due to distortion
  # newRADec = as.numeric(Rwcs_p2s(keyvalues$CRPIX1 + delta_x, keyvalues$CRPIX2 + delta_y))
  # keyvalues$CRVAL1 = newRADec[1]
  # keyvalues$CRVAL2 = newRADec[2]

  if(delta_rot != 0){
    rotmat_keyvalues = matrix(c(keyvalues$CD1_1,
                                keyvalues$CD1_2,
                                keyvalues$CD2_1,
                                keyvalues$CD2_2
    ),
    nrow = 2,
    byrow = TRUE
    )

    delta_rot_rad = -delta_rot*pi/180

    rotmat_delta = matrix(c(cos(delta_rot_rad),
                            -sin(delta_rot_rad),
                            sin(delta_rot_rad),
                            cos(delta_rot_rad)
    ),
    nrow = 2,
    byrow = TRUE
    )

    rotmat =  rotmat_keyvalues %*% rotmat_delta

    keyvalues$CD1_1 = rotmat[1,1]
    keyvalues$CD1_2 = rotmat[1,2]
    keyvalues$CD2_1 = rotmat[2,1]
    keyvalues$CD2_2 = rotmat[2,2]
  }

  if(recen){
    temp_cen = centre(keyvalues)
    keyvalues$CRVAL1 = as.numeric(temp_cen[1])
    keyvalues$CRVAL2 = as.numeric(temp_cen[2])
    keyvalues$CRPIX1 = dim(keyvalues)[1]/2 + 0.5 #the exact image centre in FITS format
    keyvalues$CRPIX2 = dim(keyvalues)[2]/2 + 0.5 #the exact image centre in FITS format
  }

  if(keylist){
    return(keyvalues)
  }else{
    input$keyvalues = keyvalues
    input = Rfits_check_image(input)

    return(input)
  }
}

.map.tran = function(x=0, y=0, delta_x=0, delta_y=0, delta_rot=0, xcen_rot=0, ycen_rot=0){
  if(delta_rot != 0){
    x_mod = (x - xcen_rot)*cos(delta_rot*pi/180) + (y - ycen_rot)*sin(delta_rot*pi/180)
    y_mod = -(x - xcen_rot)*sin(delta_rot*pi/180) + (y - ycen_rot)*cos(delta_rot*pi/180)
    x_mod = x_mod + xcen_rot + delta_x
    y_mod = y_mod + ycen_rot + delta_y
  }else{
    x_mod = x + delta_x
    y_mod = y + delta_y
  }
  list(x = x_mod, y = y_mod)
}

.cost_fn_image = function(par, image_ref, image_pre_fix, direction='backward',
                          pix_cost_use=NULL, shift_int=TRUE, return='cost', WCS_match=TRUE, cores=1L){

  if(shift_int){
    WCS_match = TRUE
  }

  if(WCS_match){
    if(inherits(image_ref, 'Rfits_image')){
      image_ref = image_ref$imDat
    }

    if(inherits(image_pre_fix, 'Rfits_image')){
      image_pre_fix = image_pre_fix$imDat
    }
  }

  if(shift_int){
    cost = .mat_diff_sum(image_ref, image_pre_fix, 1, par[1], par[2]) #1 is the unused scale param
    message(par[1],' ',par[2],' ',cost)
    return(cost)
  }else{
    if(WCS_match){
      if(length(par) == 2){
        image_post_fix = propaneTran(image_pre_fix, delta_x=par[1], delta_y=par[2], padNA=FALSE, shift_int=FALSE)
      }else if (length(par)==3){
        image_post_fix = propaneTran(image_pre_fix, delta_x=par[1], delta_y=par[2], delta_rot=par[3], padNA=FALSE, shift_int=FALSE)
      }
    }else{
      if(length(par) == 2){
        image_pre_fix = propaneWCSmod(image_pre_fix, delta_x=par[1], delta_y=par[2], recen=FALSE)
      }else{
        image_pre_fix = propaneWCSmod(image_pre_fix, delta_x=par[1], delta_y=par[2], delta_rot=par[3], recen=FALSE)
      }
      image_post_fix = propaneWarp(image_in = image_pre_fix,
                                   keyvalues_out = image_ref$keyvalues,
                                   direction = 'backward', cores=cores)$imDat
      image_ref = image_ref$imDat
    }

    if(return=='image_post_fix'){
      return(image_post_fix)
    }

    costmat = image_ref - image_post_fix

    #frac_good = length(which(!is.na(image_post_fix))) / prod(dim(image_post_fix))
    if(is.null(pix_cost_use)){
      #cost = sum(asinh((image_ref * image_post_fix)/scale/frac_good), na.rm=TRUE)
      #cost = sum(((image_ref - image_post_fix))^2, na.rm=TRUE)
      cost = sum(costmat^2, na.rm=TRUE)
    }else{
      #cost = sum(asinh((image_ref[pix_cost_use] * image_post_fix[pix_cost_use])/scale/frac_good), na.rm=TRUE)
      #cost = sum(((image_ref[pix_cost_use] - image_post_fix[pix_cost_use]))^2, na.rm=TRUE)
      cost = sum(costmat[pix_cost_use]^2, na.rm=TRUE)
    }
    #message(par[1],' ',par[2],' ',cost)
    if(return == 'cost'){
      return(cost)
    }else if(return == 'all'){
      return(list(cost=cost, image_post_fix=image_post_fix, costmat=costmat, pix_cost_use=pix_cost_use))
    }
  }
}

.cost_fn_cat = function(par, cat_ref, cat_pre_fix, rad_max, xcen_rot, ycen_rot){

  if(length(par) == 2L){
    cat_pre_fix = as.data.frame(.map.tran(cat_pre_fix[,1], cat_pre_fix[,2], delta_x=par[1], delta_y=par[2]))
  }else if(length(par) == 3L){
    cat_pre_fix = as.data.frame(.map.tran(cat_pre_fix[,1], cat_pre_fix[,2], delta_x=par[1], delta_y=par[2],
                                          delta_rot=par[3], xcen_rot=xcen_rot, ycen_rot=ycen_rot))
  }

  temp_dist = RANN::nn2(cat_pre_fix, cat_ref, searchtype='radius', radius=rad_max, k=1)$nn.dists[,1]
  temp_dist[temp_dist > 1.3e+154] = rad_max
  return(sum(temp_dist^2, na.rm=TRUE))
}

.cost_fn_image_approx = function(par, image_ref, image_pre_fix_xysub, keyvalues_pre_fix,
                                 WCS_match=TRUE, xcen_rot, ycen_rot, cores=1L){
  #image_pre_fix_xysub is just meant to be the pixels we want to transform, rather than the full image

  flux = image_pre_fix_xysub[,3]

  if(WCS_match){
    if(length(par) == 2L){
      image_pre_fix_xysub = as.data.frame(.map.tran(image_pre_fix_xysub[,1], image_pre_fix_xysub[,2], delta_x=par[1], delta_y=par[2]))
    }else if(length(par) == 3L){
      image_pre_fix_xysub = as.data.frame(.map.tran(image_pre_fix_xysub[,1], image_pre_fix_xysub[,2], delta_x=par[1], delta_y=par[2],
                                            delta_rot=par[3], xcen_rot=xcen_rot, ycen_rot=ycen_rot))
    }
  }else{
    if(length(par) == 2){
      keyvalues_pre_fix = propaneWCSmod(keyvalues_pre_fix, delta_x=par[1], delta_y=par[2], recen=FALSE)
    }else{
      keyvalues_pre_fix = propaneWCSmod(keyvalues_pre_fix, delta_x=par[1], delta_y=par[2], delta_rot=par[3], recen=FALSE)
    }
    radecsub = Rwcs_p2s(image_pre_fix_xysub[,1:2], keyvalues=keyvalues_pre_fix, cores=cores)
    image_pre_fix_xysub = Rwcs_s2p(radecsub, keyvalues=image_ref$keyvalues, cores=cores)
  }

  #xysub_loc = as.matrix(ceiling(image_pre_fix_xysub - 0.5)) #need the -0.5 because positions are in FITS format.
  #goodsel = xysub_loc[,1] >= 2L & xysub_loc[,1] <= (dim(image_ref)[1] - 1L) & xysub_loc[,2] >= 2L & xysub_loc[,2] <= (dim(image_ref)[2] - 1L)
  #image_pre_fix_xysub = image_pre_fix_xysub[goodsel,]
  #xysub_loc = xysub_loc[goodsel,]
  #flux = flux[goodsel]

  costmat = matrix(0, dim(image_ref)[1], dim(image_ref)[2])

  if(inherits(image_ref, 'Rfits_image')){
    costmat[] = image_ref$imDat
  }else{
    costmat[] = image_ref
  }

  #image_interp = matrix(0, dim(image_ref)[1], dim(image_ref)[2])

  # #cen
  # xfrac = (1 - abs(xysub_loc[,1] - 0.5 - image_pre_fix_xysub[,1]))
  # yfrac = (1 - abs(xysub_loc[,2] - 0.5 - image_pre_fix_xysub[,2]))
  # fluxshare = flux*xfrac*yfrac
  #
  # costmat[xysub_loc] = costmat[xysub_loc] - fluxshare

  #Note below we are using positions in FITS standard. This means the centre of a pixel is integer, so pixel [1,1] goes from 0.5-1.5 in x and y.

  #3x3 loops to correctly interpolate flux locally and subtract from the reference cost image
  # for(ix in -1:1){
  #   for(jy in -1:1){
  #     xfrac = (1 - abs((xysub_loc[,1] + ix) - image_pre_fix_xysub[,1]))
  #     yfrac = (1 - abs((xysub_loc[,2] + jy) - image_pre_fix_xysub[,2]))
  #
  #     if(ix != 0L | jy != 0L){ #since [0,0] is the centre pixel, which will always have some flux share between 0-1
  #       xfrac[xfrac < 0] = 0
  #       xfrac[xfrac > 1] = 0
  #       yfrac[yfrac < 0] = 0
  #       yfrac[yfrac > 1] = 0
  #     }
  #
  #     subset = cbind(xysub_loc[,1] + ix, xysub_loc[,2] + jy)
  #     subset_dup = duplicated(subset)
  #
  #     if(any(subset_dup)){
  #       costmat[subset[!subset_dup,]] = costmat[subset[!subset_dup,]] - (flux[!subset_dup]*xfrac[!subset_dup]*yfrac[!subset_dup])
  #       for(k in which(subset_dup)){
  #         costmat[subset[k,]] = costmat[subset[k,]] - (flux[k]*xfrac[k]*yfrac[k])
  #       }
  #     }else{
  #       costmat[subset] = costmat[subset] - (flux*xfrac*yfrac)
  #     }
  #   }
  # }
  #image_interp = propaneImBlur(image_interp, smooth=0.5)
  #temp = akima::interp(image_pre_fix_xysub[,1], image_pre_fix_xysub[,2], z=flux, xo=1:dim(costmat)[1], yo=1:dim(costmat)[2], linear=FALSE)$z
  #costmat = costmat - akima::interp(image_pre_fix_xysub[,1], image_pre_fix_xysub[,2], z=flux, xo=1:dim(costmat)[1], yo=1:dim(costmat)[2], linear=FALSE)$z

  #OMG, I already had this idea once... Must have been too slow!
  #cost_image = akima::interp(image_pre_fix_xysub[,1], image_pre_fix_xysub[,2], flux, xo=1:dim(image_ref)[1], yo=1:dim(image_ref)[2], linear=FALSE)$z

  # if(dump){
  #   dumpdata <<- costmat
  # }

  #costmat = costmat - image_interp

  propaneInterp2D(x = image_pre_fix_xysub[,1],
                  y = image_pre_fix_xysub[,2],
                  z = flux,
                  image = costmat,
                  type = 'sub')

  cost = sum(costmat^2, na.rm=TRUE)

  #message(paste(par, collapse=' '),' ',cost)
  return(cost)
}

.stretch = function(data, stretch='lin'){
  if(stretch == 'lin'){
    return(data)
  }else if(stretch == 'log'){
    data[which(data > 0)] = log10(data[which(data > 0)])
    return(data)
  }else if(stretch == 'asinh'){
    return(asinh(data))
  }else if(stretch == 'atan'){
    return(atan(data))
  }else if(stretch == 'sqrt'){
    return(sqrt(data))
  }else{
    stop('stretch type must be one of lin / log / asinh / atan / sqrt')
  }
}
