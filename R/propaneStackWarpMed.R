propaneStackWarpFunc = function(
    filelist = NULL,
    dirlist = NULL,
    extlist = 1,
    pattern = NULL,
    recursive = TRUE,
    zap = NULL,
    keyvalues_out = NULL,
    imager_func = NULL,
    weights = NULL,
    probs = 0.5,
    cores = floor(detectCores()/2),
    multitype = 'fork',
    chunk = 1e3,
    doweight = TRUE,
    useCUTLO = TRUE){

  timestart = proc.time()[3]

  j = NULL

  assertList(keyvalues_out, null.ok = TRUE)
  assertIntegerish(cores, len=1)
  assertIntegerish(chunk, len=1)

  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop('The Rfits package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  if(is.null(imager_func)){
    imager_func = imager::average
  }

  if(multitype=='fork'){
    registerDoParallel(cores=cores)
  }else if(multitype=='cluster'){
    registerDoParallel(cl=cores)
  }

  image_list = Rfits_make_list(filelist = filelist,
                               dirlist = dirlist,
                               extlist = extlist,
                               pattern = pattern,
                               recursive = recursive,
                               pointer = TRUE,
                               cores = cores,
                               zap = zap)

  if(!is.null(keyvalues_out)){
    which_overlap = which(foreach(i = 1:length(image_list), .combine='c')%dopar%{
      Rwcs_overlap(image_list[[i]]$keyvalues, keyvalues_ref = keyvalues_out, buffer=0)
    })
    image_list = image_list[which_overlap]
  }else{
    which_overlap = 1:length(image_list)
    temp_overlap = 1:length(image_list)
  }

  if(!is.null(keyvalues_out)){
    if(isTRUE(keyvalues_out$ZIMAGE)){
      NAXIS1 = keyvalues_out$ZNAXIS1
      NAXIS2 = keyvalues_out$ZNAXIS2
    }else{
      NAXIS1 = keyvalues_out$NAXIS1
      NAXIS2 = keyvalues_out$NAXIS2
    }
  }else{
    NAXIS1 = dim(image_list[[1]])[1]
    NAXIS2 = dim(image_list[[1]])[2]
  }

  if(NAXIS1 %% chunk == 1L){
    stop('Chunk leaves single pixel sub-array, change chunk size!')
  }

  if(NAXIS2 %% chunk == 1L){
    stop('Chunk leaves single pixel sub-array, change chunk size!')
  }

  stack_grid = expand.grid(seq(1L,NAXIS1,by=chunk), seq(1L,NAXIS2,by=chunk))
  stack_grid[,3] = stack_grid[,1] + chunk
  stack_grid[stack_grid[,3] > NAXIS1,3] = NAXIS1
  stack_grid[,4] = stack_grid[,2] + chunk
  stack_grid[stack_grid[,4] > NAXIS2,4] = NAXIS2

  stack_temp = foreach(i = 1:dim(stack_grid)[1])%dopar%{
    #for(i in 1:dim(stack_grid)[1]){
    message('Stacking sub region ',i,' of ',dim(stack_grid)[1])

    xsub = as.integer(stack_grid[i, c(1,3)])
    ysub = as.integer(stack_grid[i, c(2,4)])

    if(!is.null(keyvalues_out)){
      keyvalues_sub = Rwcs_keyvalues_sub(keyvalues_out, xsub=xsub, ysub=ysub)

      temp_overlap = which(foreach(j = 1:length(image_list), .combine = 'c')%dopar%{
        Rwcs_overlap(image_list[[j]]$keyvalues, keyvalues_sub, buffer=0)
      })
    }

    if(length(temp_overlap) == 0L){
      return(
        list(
          image = matrix(NA_real_, diff(range(xsub)) + 1L, diff(range(ysub)) + 1L),
          weight = matrix(0L, diff(range(xsub)) + 1L, diff(range(ysub)) + 1L)
        )
      )
    }

    if(!is.null(weights)){
      weights_temp = weights[temp_overlap]
      if(ifelse(is.character(imager_func), tolower(imager_func) == 'waverage', FALSE)){
        weights_temp = weights_temp/sum(weights_temp)
      }
    }

    image_list_cut = foreach(j = temp_overlap)%do%{ #not sure why this won't work in dopar... Rfits race conditions?
      if(!is.null(image_list[[j]]$keyvalues$XCUTLO) & useCUTLO){
        XCUTLO = image_list[[j]]$keyvalues$XCUTLO
      }else{
        XCUTLO = 1L #implictly assumed to be aligned images if XCUTLO is missing
      }

      if(!is.null(image_list[[j]]$keyvalues$YCUTLO) & useCUTLO){
        YCUTLO = image_list[[j]]$keyvalues$YCUTLO
      }else{
        YCUTLO = 1L #implictly assumed to be aligned images if YCUTLO is missing
      }

      xrange = c(1,(diff(range(xsub)) + 1L)) + (xsub[1] - XCUTLO)
      yrange = c(1,(diff(range(ysub)) + 1L)) + (ysub[1] - YCUTLO)
      return(imager::as.cimg(image_list[[j]][xrange, yrange]$imDat))
    }

    if(is.function(imager_func)){
      if('w' %in% names(formals(imager_func))){
        if(is.null(weights)){
          stop('weights must be specified for imager wsum')
        }
        image_stack = as.matrix(imager_func(image_list_cut, w=weights_temp, na.rm=TRUE)) #only relevant for wsum
      }else if('prob' %in% names(formals(imager_func))){
        image_stack = as.matrix(imager_func(image_list_cut, prob=probs, na.rm=TRUE))
      }else if('increasing' %in% names(formals(imager_func))){
        image_stack = imager_func(image_list, increasing=TRUE)
      }else{
        if('na.rm' %in% names(formals(imager_func))){
          image_stack = as.matrix(imager_func(image_list_cut, na.rm=TRUE))
        }else{
          image_stack = as.matrix(imager_func(image_list_cut))
        }
      }
    }else{
      if(tolower(imager_func) == 'invar'){
        image_stack = 1/as.matrix(imager::parvar(image_list_cut, na.rm=TRUE))
      }else if(tolower(imager_func) == 'quantile' | tolower(imager_func) == 'quan' ){
        if(!requireNamespace("ProUtils", quietly = TRUE)){
          stop('The ProUtils package is needed to compute quantiles, install from asgr/ProUtils')
        }

        im_dim = dim(image_list_cut[[1]])
        temp_mat = matrix(as.numeric(unlist(image_list_cut)), nrow=length(image_list_cut), byrow = TRUE)

        if(is.null(weights)){
          temp_weight = NULL
        }else{
          temp_weight = matrix(as.numeric(unlist(weights)), nrow=length(image_list_cut), byrow = TRUE)
        }

        temp_out = ProUtils::quan_wt_mat_col(temp_mat, probs=probs, wt=temp_weight)

        image_stack = list()

        if(length(probs) == 1){
          image_stack = matrix(temp_out, im_dim[1], im_dim[2])
        }else{
          for(k in 1:length(probs)){
            image_stack = c(image_stack, list(matrix(temp_out[,k], im_dim[1], im_dim[2])))
          }
        }
        # old code (slow!)
        # image_stack = list()
        # for(k in 1:length(probs)){
        #   temp_out = apply(temp_mat, MARGIN=2, FUN=quantile, probs=probs[k], na.rm=TRUE)
        #   image_stack = c(image_stack, list(matrix(temp_out, im_dim[1], im_dim[2])))
        # }
      }else if(tolower(imager_func) == 'waverage'){
        image_stack = as.matrix(imager::wsum(image_list_cut, w=weights_temp, na.rm=TRUE))
      }
    }

    if(doweight | ifelse(is.character(imager_func), tolower(imager_func) == 'invar', FALSE)){
      weight_list = foreach(j = 1:length(image_list_cut))%do%{
        return(imager::as.cimg(!is.na(image_list_cut[[j]])))
      }

      weight_stack = as.matrix(imager::add(weight_list))

      if(ifelse(is.character(imager_func), tolower(imager_func) == 'invar', FALSE)){
        image_stack = image_stack*(weight_stack - 1)
        image_stack[weight_stack < 2] = NA
      }
    }else{
      weight_stack = NULL
    }

    return(
      list(image = image_stack, weight = weight_stack)
    )
  }

  image_stack = matrix(0, NAXIS1, NAXIS2)
  if(doweight){
    weight_stack = matrix(0L, NAXIS1, NAXIS2)
  }else{
    weight_stack = NULL
  }

  for(i in 1:dim(stack_grid)[1]){
    image_stack[stack_grid[i,1]:stack_grid[i,3], stack_grid[i,2]:stack_grid[i,4]] = stack_temp[[i]]$image
    if(doweight){
      weight_stack[stack_grid[i,1]:stack_grid[i,3], stack_grid[i,2]:stack_grid[i,4]] = stack_temp[[i]]$weight
    }
  }

  if(!is.null(keyvalues_out)){
    keyvalues_out$EXTNAME = 'image'
    keyvalues_out$MAGZERO = image_list[[1]]$keyvalues$MAGZERO
    keyvalues_out$R_VER = R.version$version.string
    keyvalues_out$PANE_VER = as.character(packageVersion('ProPane'))
    keyvalues_out$RWCS_VER = as.character(packageVersion('Rwcs'))

    image_stack = Rfits_create_image(image_stack,
                                     keyvalues=keyvalues_out,
                                     keypass=FALSE,
                                     history='Stacked with Rwcs_stack')

    if(doweight){
      keyvalues_out$EXTNAME = 'weight'
      keyvalues_out$MAGZERO = NULL
      weight_stack = Rfits_create_image(weight_stack,
                                        keyvalues=keyvalues_out,
                                        keypass=FALSE)
    }
  }

  time_taken = proc.time()[3] - timestart
  message('Time taken: ',signif(time_taken,4),' seconds')

  output = list(
    image = image_stack,
    weight = weight_stack,
    which_overlap = which_overlap,
    time = time_taken,
    Nim = length(which_overlap)
  )

  class(output) = "ProPane"
  return(invisible(output))
}

propaneStackWarpMed = function(
      filelist = NULL,
      dirlist = NULL,
      extlist = 1,
      pattern = NULL,
      recursive = TRUE,
      zap = NULL,
      keyvalues_out = NULL,
      cores = floor(detectCores()/2),
      multitype = 'fork',
      chunk = 1e3,
      doweight = TRUE,
      useCUTLO = TRUE){

  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  output = propaneStackWarpFunc(
    filelist = filelist,
    dirlist = dirlist,
    extlist = extlist,
    pattern = pattern,
    recursive = recursive,
    zap = zap,
    imager_func = imager::parmed,
    keyvalues_out = keyvalues_out,
    probs = 0.5,
    cores = cores,
    multitype = multitype,
    chunk = chunk,
    doweight = doweight,
    useCUTLO = useCUTLO
  )

  return(output)
}
