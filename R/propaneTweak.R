propaneTweak = function(image_ref, image_pre_fix, delta_max=c(3,0), quan_cut=0.99, Nmeta=3,
                      WCS_match=TRUE, cores=4, shift_int=TRUE, return_image=TRUE, direction='backward',
                      final_centre=TRUE, cutcheck=FALSE, verbose=TRUE){

  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
  }

  timestart = proc.time()[3]

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

  dim_orig = dim(image_ref)

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
  }

  if(inherits(image_pre_fix, 'Rfits_image')){
    if(WCS_match){
      im_med = median(image_pre_fix$imDat, na.rm=TRUE)
      im_quan_lo = quantile(image_pre_fix$imDat, quan_cut[1], na.rm=TRUE)
    }else{
      image_ref_warp = propaneWarp(image_ref,
                                   keyvalues_out = image_pre_fix$keyvalues,
                                   direction = 'backward'
      )$imDat
      sel = which(!is.na(image_ref_warp))
      im_med = median(image_pre_fix$imDat[sel], na.rm=TRUE)
      im_quan_lo = quantile(image_pre_fix$imDat[sel], quan_cut[1], na.rm=TRUE)
    }
    if(length(quan_cut)==2){
      image_pre_fix$imDat[image_pre_fix$imDat > quantile(image_pre_fix$imDat, quan_cut[2], na.rm=TRUE)] = NA
    }
    image_pre_fix$imDat = image_pre_fix$imDat - im_med
    image_pre_fix$imDat = image_pre_fix$imDat / im_quan_lo
    image_pre_fix$imDat[image_pre_fix$imDat < 1] = NA
    if(cutcheck){
      plot(image_pre_fix)
      legend('topleft', 'image_pre_fix')
    }
  }else{
    if(length(quan_cut)==2){
      image_pre_fix[image_pre_fix > quantile(image_pre_fix, quan_cut[2], na.rm=TRUE)] = NA
    }
    image_pre_fix = image_pre_fix - median(image_pre_fix, na.rm=TRUE)
    image_pre_fix = image_pre_fix / quantile(image_pre_fix, quan_cut[1], na.rm=TRUE)
    image_pre_fix[image_pre_fix < 1] = NA
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
    if(cutcheck){
      magimage(image_ref)
      legend('topleft', 'image_ref')
      return(NULL)
    }
  }


  scale = 1 # not really used anymore
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
        cost = .mat_diff_sum(image_ref, image_pre_fix, scale, grid_search[i,1], grid_search[i,2])
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
            cost = .mat_diff_sum(image_ref, image_pre_fix, scale, grid_search[i,1], grid_search[i,2])
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

    optim_out = optim(par = par,
                      fn = .cost_fn,
                      method = "L-BFGS-B",
                      lower = lower,
                      upper = upper,
                      image_ref = image_ref,
                      image_pre_fix = image_pre_fix,
                      scale = scale,
                      direction = direction,
                      pix_cost_use = NULL,
                      shift_int = FALSE,
                      WCS_match = WCS_match
    )
  }

  if(return_image){
    if(verbose){
      message('Creating output image_post_fix')
    }

    if(inherits(image_pre_fix_orig, 'Rfits_image')){
      if(WCS_match){
        image_post_fix = image_pre_fix_orig

        if(optim_out$par[1] != 0 | optim_out$par[2] != 0){
          image_post_fix$imDat = .cost_fn(par = optim_out$par,
                                    image_ref = image_ref,
                                    image_pre_fix = image_pre_fix_orig$imDat,
                                    scale = scale,
                                    direction = direction,
                                    return = 'image_post_fix',
                                    shift_int = FALSE,
                                    WCS_match = TRUE)
        }
      }else{
        if(length(optim_out$par) == 2){
          image_post_fix = propaneWCSmod(image_pre_fix_orig, delta_x=optim_out$par[1], delta_y=optim_out$par[2])
        }else{
          image_post_fix = propaneWCSmod(image_pre_fix_orig, delta_x=optim_out$par[1], delta_y=optim_out$par[2], delta_rot=optim_out$par[3])
        }
      }
    }else{
      if(optim_out$par[1] != 0 | optim_out$par[2] != 0){
        image_post_fix = .cost_fn(par = optim_out$par,
                                       image_ref = image_ref,
                                       image_pre_fix = image_pre_fix_orig,
                                       scale = scale,
                                       direction = direction,
                                       return = 'image_post_fix',
                                       shift_int = FALSE,
                                       WCS_match = TRUE)
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

propaneTran = function(image, delta_x = 0, delta_y = 0, delta_rot = 0, xcen_rot = dim(image)[1]/2,
                       ycen_rot = dim(image)[1]/2, direction='backward', padNA=TRUE, shift_int=TRUE){
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

propaneWCSmod = function(input, delta_x = 0, delta_y = 0, delta_rot = 0){

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

  header = Rfits_keyvalues_to_header(keyvalues)
  hdr = Rfits_keyvalues_to_hdr(keyvalues)
  raw = Rfits_header_to_raw(header)

  if(keylist){
    return(list(keyvalues = keyvalues,
                header = header,
                hdr = hdr,
                raw = raw)
    )
  }else{
    input$keyvalues = keyvalues
    input$header = header
    input$hdr = hdr
    input$raw = raw

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

.cost_fn = function(par, image_ref, image_pre_fix, scale=1, direction='backward', pix_cost_use=NULL,
                    shift_int=TRUE, return='cost', WCS_match=TRUE){

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
    cost = .mat_diff_sum(image_ref, image_pre_fix, scale, par[1], par[2])
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
        image_pre_fix = propaneWCSmod(image_pre_fix, delta_x=par[1], delta_y=par[2])
      }else{
        image_pre_fix = propaneWCSmod(image_pre_fix, delta_x=par[1], delta_y=par[2], delta_rot=par[3])
      }
      image_post_fix = propaneWarp(image_in = image_pre_fix,
                                   keyvalues_out = image_ref$keyvalues,
                                   direction = 'backward')$imDat
      image_ref = image_ref$imDat
    }
    if(return=='image_post_fix'){
      return(image_post_fix)
    }
    #frac_good = length(which(!is.na(image_post_fix))) / prod(dim(image_post_fix))
    if(is.null(pix_cost_use)){
      #cost = sum(asinh((image_ref * image_post_fix)/scale/frac_good), na.rm=TRUE)
      cost = sum(((image_ref - image_post_fix)/scale)^2, na.rm=TRUE)
    }else{
      #cost = sum(asinh((image_ref[pix_cost_use] * image_post_fix[pix_cost_use])/scale/frac_good), na.rm=TRUE)
      cost = sum(((image_ref - image_post_fix[pix_cost_use])/scale)^2, na.rm=TRUE)
    }
    #message(par[1],' ',par[2],' ',cost)
    if(return=='cost'){
      return(cost)
    }
    if(return=='all'){
      return(list(cost=cost, image_post_fix=image_post_fix))
    }
  }
}
